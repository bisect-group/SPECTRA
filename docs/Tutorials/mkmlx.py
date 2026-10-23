#!/usr/bin/env python3
"""Build a MATLAB live script (.mlx) from a plain-text source file.

The .mlx format is an OPC zip holding WordprocessingML, so it cannot be
written or diffed by hand. This builds one from a readable text source, which
means tutorials can be written and reviewed as text and still ship as live
scripts.

Source format, one block per paragraph, blocks separated by blank lines::

    # Title                  ->  title style, once, at the top
    ## Heading               ->  heading style
    ### Subheading           ->  heading2 style
    ```                      ->  code style; everything up to the closing
    code goes here               fence is kept verbatim, including blank
    ```                          lines and indentation
    Anything else            ->  text style (a paragraph; wrapped lines in
                                 the source are joined into one)

Usage::

    python3 mkmlx.py Tutorial_consistency.txt            # -> Tutorial_consistency.mlx
    python3 mkmlx.py src/*.txt -o .

Note:
    The package parts other than the document body are fixed boilerplate
    copied from a live script MATLAB itself wrote. They are not guessed:
    MATLAB rejects the file outright ("Can't open file") if the relationship
    types or namespaces are even slightly wrong.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys
import zipfile

CONTENT_TYPES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>'
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
    '<Default ContentType="application/vnd.openxmlformats-package.relationships+xml"'
    ' Extension="rels"/>'
    '<Default ContentType="application/vnd.mathworks.matlab.code.document+xml"'
    ' Extension="xml"/>'
    '<Override ContentType="text/xml" PartName="/matlab/output.xml"/>'
    '<Override ContentType="application/vnd.openxmlformats-package.core-properties+xml"'
    ' PartName="/metadata/coreProperties.xml"/>'
    '<Override ContentType="application/vnd.mathworks.package.coreProperties+xml"'
    ' PartName="/metadata/mwcoreProperties.xml"/>'
    '<Override ContentType="application/vnd.mathworks.package.corePropertiesExtension+xml"'
    ' PartName="/metadata/mwcorePropertiesExtension.xml"/>'
    "</Types>"
)

RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Target="matlab/document.xml"'
    ' Type="http://schemas.mathworks.com/matlab/code/2013/relationships/document"/>'
    '<Relationship Id="rId2" Target="matlab/output.xml"'
    ' Type="http://schemas.mathworks.com/matlab/code/2013/relationships/output"/>'
    '<Relationship Id="rId3" Target="metadata/mwcoreProperties.xml"'
    ' Type="http://schemas.mathworks.com/package/2012/relationships/coreProperties"/>'
    '<Relationship Id="rId4" Target="metadata/mwcorePropertiesExtension.xml"'
    ' Type="http://schemas.mathworks.com/package/2014/relationships/corePropertiesExtension"/>'
    '<Relationship Id="rId5" Target="metadata/coreProperties.xml"'
    ' Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties"/>'
    "</Relationships>"
)

OUTPUT = (
    '<?xml version="1.0" encoding="UTF-8"?>'
    "<embeddedOutputs><metaData>"
    "<evaluationState>manual</evaluationState>"
    "<layoutState>document</layoutState>"
    "<outputStatus>ready</outputStatus>"
    "</metaData>"
    '<outputArray type="array"/><regionArray type="array"/>'
    "</embeddedOutputs>"
)

CORE_PROPERTIES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>'
    '<cp:coreProperties'
    ' xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"'
    ' xmlns:dc="http://purl.org/dc/elements/1.1/"'
    ' xmlns:dcmitype="http://purl.org/dc/dcmitype/"'
    ' xmlns:dcterms="http://purl.org/dc/terms/"'
    ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    '<dcterms:created xsi:type="dcterms:W3CDTF">2016-07-14T14:56:22Z</dcterms:created>'
    '<dcterms:modified xsi:type="dcterms:W3CDTF">2017-05-01T14:39:40Z</dcterms:modified>'
    "</cp:coreProperties>"
)

MW_CORE_PROPERTIES = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>'
    '<mwcoreProperties xmlns="http://schemas.mathworks.com/package/2012/coreProperties">'
    "<contentType>application/vnd.mathworks.matlab.code</contentType>"
    "<contentTypeFriendlyName>MATLAB Code</contentTypeFriendlyName>"
    "<matlabRelease>R2017a</matlabRelease>"
    "</mwcoreProperties>"
)

MW_CORE_PROPERTIES_EXT = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>'
    '<mwcoreProperties'
    ' xmlns="http://schemas.mathworks.com/package/2014/corePropertiesExtension">'
    "<matlabVersion>9.0.0.383818</matlabVersion>"
    "</mwcoreProperties>"
)

FENCE = re.compile(r"^```\s*$")


def escape(text: str) -> str:
    """Escape text for an XML text node."""
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def parse(source: str) -> list[tuple[str, str]]:
    """Turn the text source into a list of (style, text) paragraphs."""
    blocks: list[tuple[str, str]] = []
    lines = source.splitlines()
    i = 0
    paragraph: list[str] = []

    def flush() -> None:
        if paragraph:
            blocks.append(("text", " ".join(paragraph).strip()))
            paragraph.clear()

    while i < len(lines):
        line = lines[i]

        if FENCE.match(line):
            flush()
            i += 1
            code: list[str] = []
            while i < len(lines) and not FENCE.match(lines[i]):
                code.append(lines[i])
                i += 1
            if i >= len(lines):
                raise SystemExit("unterminated code fence")
            i += 1
            while code and not code[0].strip():
                code.pop(0)
            while code and not code[-1].strip():
                code.pop()
            blocks.append(("code", "\n".join(code)))
            continue

        if not line.strip():
            flush()
        elif line.startswith("### "):
            flush()
            blocks.append(("heading2", line[4:].strip()))
        elif line.startswith("## "):
            flush()
            blocks.append(("heading", line[3:].strip()))
        elif line.startswith("# "):
            flush()
            blocks.append(("title", line[2:].strip()))
        else:
            paragraph.append(line.strip())
        i += 1

    flush()
    return blocks


def document(blocks: list[tuple[str, str]]) -> str:
    """Render paragraphs as the WordprocessingML document body."""
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        "<w:body>"
    ]
    for style, text in blocks:
        # Code keeps its newlines and indentation, so it goes in a CDATA
        # section; prose is escaped normally.
        if style == "code":
            if "]]>" in text:
                raise SystemExit("code block contains ]]>, which CDATA cannot hold")
            body = f"<![CDATA[{text}]]>"
        else:
            body = escape(text)
        parts.append(
            f'<w:p><w:pPr><w:pStyle w:val="{style}"/></w:pPr>'
            f"<w:r><w:t>{body}</w:t></w:r></w:p>"
        )
    parts.append("</w:body></w:document>")
    return "".join(parts)


def build(source_path: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    blocks = parse(source_path.read_text(encoding="utf-8"))
    if not blocks:
        raise SystemExit(f"{source_path}: nothing to write")

    target = out_dir / (source_path.stem + ".mlx")
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as package:
        package.writestr("[Content_Types].xml", CONTENT_TYPES)
        package.writestr("_rels/.rels", RELS)
        package.writestr("matlab/document.xml", document(blocks))
        package.writestr("matlab/output.xml", OUTPUT)
        package.writestr("metadata/coreProperties.xml", CORE_PROPERTIES)
        package.writestr("metadata/mwcoreProperties.xml", MW_CORE_PROPERTIES)
        package.writestr("metadata/mwcorePropertiesExtension.xml", MW_CORE_PROPERTIES_EXT)

    counts = {}
    for style, _ in blocks:
        counts[style] = counts.get(style, 0) + 1
    summary = ", ".join(f"{n} {style}" for style, n in sorted(counts.items()))
    print(f"{target.name}: {summary}")
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sources", nargs="+", type=pathlib.Path)
    parser.add_argument("-o", "--out-dir", type=pathlib.Path, default=None,
                        help="where to write the .mlx (default: next to the source)")
    args = parser.parse_args()

    for source in args.sources:
        if not source.is_file():
            print(f"{source}: not a file", file=sys.stderr)
            return 1
        build(source, args.out_dir or source.parent)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
