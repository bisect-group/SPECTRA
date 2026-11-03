"""Attach the toy-network panel to the four-arm CC figure.

The main-text figure is not the bare bar chart: each of the three runtime
panels has a toy network pasted on its left and a caption strip on top. That
is what merge_runtime_images.ipynb does, and its cell 1 builds exactly this
panel - consistency_bg/Slide2.PNG + SPECTRA_CC_runtime.png ->
SPECTRA_CC_runtime_merged.png, captioned "Efficiency comparison for
consistency check (CC)".

This script reuses those two functions rather than reimplementing them: it
reads the notebook, executes only the cell that defines merge_side_by_side
and image_merge, and calls them with the new four-arm chart in place of the
old three-arm one. Nothing about the layout, padding, caption font or height
changes, so the regenerated panel drops into the figure beside the untouched
SPECTRA_ME and SPECTRA_CC_ME panels.
"""

import contextlib
import json
import os
import shutil
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.join(HERE, "merge_runtime_images.ipynb")


@contextlib.contextmanager
def dejavu_on_path():
    """Make ImageFont.truetype('DejaVuSans-Bold.ttf') resolve.

    image_merge asks PIL for that font by bare filename. PIL looks in the
    working directory and a few system font directories, none of which hold
    DejaVu on macOS, so the call raises "cannot open resource" - on Linux,
    where the notebook was written, it is usually installed system-wide.

    Rather than edit the notebook's function, drop a copy of the font - the
    one matplotlib ships, so no download - into a temporary directory and run
    from there. The output paths are absolute, so nothing lands in the temp
    directory but the font, and it goes away with it.
    """
    import matplotlib

    src = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf")
    if not os.path.exists(src):
        raise FileNotFoundError(f"matplotlib has no DejaVuSans-Bold.ttf at {src}")
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as td:
        shutil.copy(src, os.path.join(td, "DejaVuSans-Bold.ttf"))
        os.chdir(td)
        try:
            yield
        finally:
            os.chdir(cwd)


def load_notebook_functions():
    """Execute the notebook cell that defines the two merge helpers."""
    nb = json.load(open(NB))
    ns = {}
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        # The defining cell is the one carrying both functions; the calling
        # cells are skipped so importing this does not regenerate the other
        # two panels as a side effect.
        if "def merge_side_by_side" in src and "def image_merge" in src:
            exec(compile(src, NB, "exec"), ns)
            return ns["merge_side_by_side"], ns["image_merge"]
    raise RuntimeError(f"no cell in {NB} defines both merge helpers")


def main():
    merge_side_by_side, image_merge = load_notebook_functions()

    # Same left panel and same caption as the published CC figure; only the
    # right-hand chart differs.
    im1_path = os.path.join(HERE, "consistency_bg", "Slide2.PNG")
    im2_path = os.path.join(HERE, "SPECTRA_CC_FASTERCC_runtime.png")
    out_path = os.path.join(HERE, "SPECTRA_CC_FASTERCC_runtime_merged.png")

    merged_path = merge_side_by_side(
        im1_path, im2_path, out_path=out_path, match_to="right", padding=10
    )
    print("Merged image saved to:", merged_path)

    with dejavu_on_path():
        image_merge(
            out_path,
            out_path,
            caption_text="Efficiency comparison for consistency check (CC)",
            caption_font_size=50,
            caption_height=100,
        )


if __name__ == "__main__":
    main()
