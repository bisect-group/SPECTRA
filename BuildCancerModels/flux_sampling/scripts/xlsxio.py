"""Minimal reader/writer for the .xlsx files in this project.

openpyxl is not installed in the cluster scipy-stack, and the archived FVA
sheets use inline strings rather than a shared-string table, so a small
purpose-built parser is simpler than adding a dependency.
"""
import re
import zipfile
from xml.sax.saxutils import escape

# Attributes appear in different orders depending on the writer, so the tag is
# matched whole and its attributes are pulled out individually rather than
# assuming a fixed layout.
_CELL = re.compile(rb'<c\b([^>]*?)(?:/>|>(.*?)</c>)', re.S)
_REF = re.compile(rb'\br="([A-Z]+)(\d+)"')
_TYPE = re.compile(rb'\bt="([^"]+)"')
_V = re.compile(rb'<v[^>]*>(.*?)</v>', re.S)
_T = re.compile(rb'<t[^>]*>(.*?)</t>', re.S)


def _col_to_idx(col):
    n = 0
    for ch in col:
        n = n * 26 + (ch - 64)
    return n - 1


def _unescape(b):
    s = b.decode('utf-8', 'ignore')
    return (s.replace('&lt;', '<').replace('&gt;', '>')
             .replace('&quot;', '"').replace('&apos;', "'")
             .replace('&amp;', '&'))


def read_sheets(path):
    """Return {sheet_name: [[cell, ...], ...]} for every worksheet."""
    z = zipfile.ZipFile(path)
    names = z.namelist()

    shared = []
    if 'xl/sharedStrings.xml' in names:
        shared = [_unescape(m.group(1))
                  for m in _T.finditer(z.read('xl/sharedStrings.xml'))]

    wb = z.read('xl/workbook.xml').decode('utf-8', 'ignore')
    rels = z.read('xl/_rels/workbook.xml.rels').decode('utf-8', 'ignore')
    relmap = {}
    for tag in re.findall(r'<Relationship\b[^>]*/?>', rels):
        rid = re.search(r'\bId="([^"]+)"', tag)
        tgt = re.search(r'\bTarget="([^"]+)"', tag)
        if rid and tgt:
            relmap[rid.group(1)] = tgt.group(1)
    sheets = []
    for tag in re.findall(r'<sheet\b[^>]*/?>', wb):
        nm = re.search(r'\bname="([^"]*)"', tag)
        rid = re.search(r'\br:id="([^"]+)"', tag)
        if nm and rid:
            sheets.append((nm.group(1), rid.group(1)))

    out = {}
    for name, rid in sheets:
        target = relmap[rid].lstrip('/')
        if not target.startswith('xl/'):
            target = 'xl/' + target
        if target not in names:
            continue
        xml = z.read(target)
        rows = {}
        for m in _CELL.finditer(xml):
            attrs, body = m.group(1), m.group(2) or b''
            ref = _REF.search(attrs)
            if not ref:
                continue
            col, rownum = ref.group(1), int(ref.group(2))
            tm = _TYPE.search(attrs)
            ctype = tm.group(1).decode() if tm else 'n'
            if ctype == 'inlineStr':
                tm = _T.search(body)
                val = _unescape(tm.group(1)) if tm else ''
            elif ctype == 's':
                vm = _V.search(body)
                val = shared[int(vm.group(1))] if vm else ''
            elif ctype == 'str':
                vm = _V.search(body)
                val = _unescape(vm.group(1)) if vm else ''
            else:
                vm = _V.search(body)
                if not vm or not vm.group(1):
                    val = None
                else:
                    try:
                        val = float(vm.group(1))
                    except ValueError:
                        val = _unescape(vm.group(1))
            rows.setdefault(rownum, {})[_col_to_idx(col)] = val
        if not rows:
            out[name] = []
            continue
        maxcol = max(max(r) for r in rows.values()) + 1
        out[name] = [[rows.get(i, {}).get(j) for j in range(maxcol)]
                     for i in range(1, max(rows) + 1)]
    return out


def write_sheets(path, sheets):
    """sheets: {name: [[cell, ...], ...]}. Everything is written as inline text
    or a number, which is all the downstream consumers here need."""
    def cell_xml(ref, v):
        if v is None or v == '':
            return ''
        if isinstance(v, bool):
            v = int(v)
        if isinstance(v, (int, float)):
            return f'<c r="{ref}" t="n"><v>{v!r}</v></c>' if isinstance(v, float) \
                   else f'<c r="{ref}" t="n"><v>{v}</v></c>'
        return f'<c r="{ref}" t="inlineStr"><is><t xml:space="preserve">{escape(str(v))}</t></is></c>'

    def colname(i):
        s = ''
        i += 1
        while i:
            i, r = divmod(i - 1, 26)
            s = chr(65 + r) + s
        return s

    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml',
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            + ''.join(f'<Override PartName="/xl/worksheets/sheet{i+1}.xml" '
                      'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                      for i in range(len(sheets)))
            + '</Types>')
        z.writestr('_rels/.rels',
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '</Relationships>')
        z.writestr('xl/workbook.xml',
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
            + ''.join(f'<sheet name="{escape(n)[:31]}" sheetId="{i+1}" r:id="rId{i+1}"/>'
                      for i, n in enumerate(sheets))
            + '</sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels',
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + ''.join(f'<Relationship Id="rId{i+1}" '
                      'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                      f'Target="worksheets/sheet{i+1}.xml"/>' for i in range(len(sheets)))
            + '</Relationships>')
        for i, (name, rows) in enumerate(sheets.items()):
            body = []
            for ri, row in enumerate(rows, start=1):
                cells = ''.join(cell_xml(f'{colname(ci)}{ri}', v) for ci, v in enumerate(row))
                body.append(f'<row r="{ri}">{cells}</row>')
            z.writestr(f'xl/worksheets/sheet{i+1}.xml',
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                '<sheetData>' + ''.join(body) + '</sheetData></worksheet>')
