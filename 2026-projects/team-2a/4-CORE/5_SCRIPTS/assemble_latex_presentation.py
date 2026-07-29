#!/usr/bin/env python3

from pathlib import Path
import sys
import subprocess

def escape_tex(s):
    """Escape characters that are special in latex."""
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }

    for old, new in replacements.items():
        s = s.replace(old, new)
    return s


def main(image_dir, output_file):
    image_dir = Path(image_dir)

    images = sorted(image_dir.glob("*.png"))

    if not images:
        raise RuntimeError(f"No PNG files found in '{image_dir}'")

    # write a beamer latex file
    with open(image_dir/output_file, "w") as f:
        f.write(r"""\documentclass{beamer}

\usepackage{graphicx}

\setbeamertemplate{navigation symbols}{}

\begin{document}

""")

        for img in images:
            title = escape_tex(img.stem)
            #subtitle = str(image_dir.resolve())
            subtitle = str(image_dir.resolve()).split("/", -1)[-1]
            subtitle = subtitle.replace("_", "\_")

            f.write(r"\begin{frame}{" + title + "}{" + subtitle + "}\n")
            f.write(r"\centering" + "\n")
            f.write(
                r"\includegraphics[width=\textwidth,height=0.9\textheight,keepaspectratio]{"
                + img.as_posix()
                + "}\n"
            )
            f.write(r"\end{frame}" + "\n\n")

        f.write(r"\end{document}")

    # assemble pdf
    #build = Path("pdf")
    build = image_dir/"pdf"
    build.mkdir(exist_ok=True)
    try:
        subprocess.run(
            [
                "pdflatex",
                "-interaction=nonstopmode",
                f"-output-directory={build}",
                image_dir/output_file
            ],
            check=True
        )
    except FileNotFoundError as e:
        print("\n[PDF Generation skipped]")
        print("Reason: pdflatex compiler is not available on this system.")
        print(f"Detail: {e}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:")
        print(f"\tpython {sys.argv[0]} image_directory presentation.tex")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
