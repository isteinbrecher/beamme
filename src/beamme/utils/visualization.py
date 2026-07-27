# The MIT License (MIT)
#
# Copyright (c) 2018-2026 BeamMe Authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Helper functions for visualizations."""

import base64 as _base64
import html as _html
import os as _os
import tempfile as _tempfile
import uuid as _uuid
from pathlib import Path as _Path
from typing import Sequence as _Sequence

import ipywidgets as _widgets
import pyvista as _pv
from IPython.display import HTML as _HTML
from IPython.display import IFrame as _IFrame
from IPython.display import display as _display

from beamme.utils.environment import is_nbsphinx as _is_nbsphinx


def _validate_plotter_sequence(
    plotters: _Sequence[_pv.Plotter],
    labels: _Sequence[str] | None,
) -> list[str]:
    """Validate a plotter sequence and return its frame labels."""
    if len(plotters) == 0:
        raise ValueError("At least one plotter is required.")

    if labels is None:
        return [f"Frame {index + 1}" for index in range(len(plotters))]

    if len(labels) != len(plotters):
        raise ValueError("The number of labels must match the number of plotters.")

    return [str(label) for label in labels]


def _sequence_viewer_html(
    static_sources: _Sequence[str],
    labels: _Sequence[str],
    *,
    interactive_sources: _Sequence[str] | None,
) -> str:
    """Create the HTML for a static and interactive plotter sequence viewer."""
    viewer_uid = f"beamme-sequence-{_uuid.uuid4().hex}"
    escaped_labels = [_html.escape(label) for label in labels]

    def create_frames(sources, frame_type):
        frames = []
        for index, source in enumerate(sources):
            hidden = "" if index == 0 else " hidden"
            escaped_source = _html.escape(source, quote=True)
            escaped_label = escaped_labels[index]
            if frame_type == "static":
                content = (
                    f'<img src="{escaped_source}" alt="{escaped_label}" '
                    'style="width:100%; height:100%; object-fit:contain;">'
                )
            else:
                content = (
                    f'<iframe src="{escaped_source}" title="{escaped_label}" '
                    'style="width:100%; height:100%; border:0;"></iframe>'
                )
            frames.append(
                f'<div class="beamme-sequence-frame"{hidden}>{content}</div>'
            )
        return "".join(frames)

    static_frames = create_frames(static_sources, "static")
    interactive_tab = ""
    interactive_panel = ""
    if interactive_sources is not None:
        interactive_frames = create_frames(interactive_sources, "interactive")
        interactive_tab = """
            <button type="button" role="tab" data-view="interactive"
                aria-selected="false">Interactive 3D</button>
        """
        interactive_panel = f"""
            <div class="beamme-sequence-panel" data-panel="interactive" hidden>
                {interactive_frames}
            </div>
        """

    labels_javascript = ", ".join(
        f'"{_html.escape(label, quote=True)}"' for label in labels
    )
    return f"""
    <div id="{viewer_uid}" class="beamme-sequence">
        <style>
            #{viewer_uid} {{
                border: 1px solid #d6d6d6;
                border-radius: 8px;
                overflow: hidden;
                background: #fff;
                color: #222;
                font-family: sans-serif;
            }}
            #{viewer_uid} .beamme-sequence-tabs {{
                display: flex;
                gap: 0.25rem;
                padding: 0.5rem 0.5rem 0;
            }}
            #{viewer_uid} button {{
                border: 1px solid #aaa;
                border-radius: 4px;
                padding: 0.35rem 0.7rem;
                background: #f4f4f4;
                cursor: pointer;
            }}
            #{viewer_uid} [role="tab"][aria-selected="true"] {{
                background: #ddd;
                font-weight: 600;
            }}
            #{viewer_uid} .beamme-sequence-panel {{
                height: 540px;
                padding: 0.5rem;
            }}
            #{viewer_uid} .beamme-sequence-frame {{
                width: 100%;
                height: 100%;
            }}
            #{viewer_uid} .beamme-sequence-controls {{
                display: flex;
                align-items: center;
                justify-content: center;
                gap: 0.75rem;
                padding: 0.5rem;
                border-top: 1px solid #d6d6d6;
            }}
            #{viewer_uid} output {{
                min-width: 12rem;
                text-align: center;
            }}
        </style>
        <div class="beamme-sequence-tabs" role="tablist">
            <button type="button" role="tab" data-view="static"
                aria-selected="true">Static frames</button>
            {interactive_tab}
        </div>
        <div class="beamme-sequence-panel" data-panel="static">
            {static_frames}
        </div>
        {interactive_panel}
        <div class="beamme-sequence-controls">
            <button type="button" data-step="-1" aria-label="Previous frame">
                &#9664; Previous
            </button>
            <output aria-live="polite"></output>
            <button type="button" data-step="1" aria-label="Next frame">
                Next &#9654;
            </button>
        </div>
        <script>
            (() => {{
                const root = document.getElementById("{viewer_uid}");
                const labels = [{labels_javascript}];
                let frameIndex = 0;

                function updateFrames() {{
                    root.querySelectorAll(".beamme-sequence-panel").forEach(
                        panel => {{
                            panel.querySelectorAll(
                                ".beamme-sequence-frame"
                            ).forEach((frame, index) => {{
                                frame.hidden = index !== frameIndex;
                            }});
                        }}
                    );
                    root.querySelector("output").textContent =
                        `${{frameIndex + 1}} / ${{labels.length}} — ` +
                        labels[frameIndex];
                }}

                root.querySelectorAll("[data-step]").forEach(button => {{
                    button.addEventListener("click", () => {{
                        const step = Number(button.dataset.step);
                        frameIndex =
                            (frameIndex + step + labels.length) % labels.length;
                        updateFrames();
                    }});
                }});

                root.querySelectorAll("[role=tab]").forEach(tab => {{
                    tab.addEventListener("click", () => {{
                        root.querySelectorAll("[role=tab]").forEach(item => {{
                            item.setAttribute(
                                "aria-selected", item === tab ? "true" : "false"
                            );
                        }});
                        root.querySelectorAll(
                            ".beamme-sequence-panel"
                        ).forEach(panel => {{
                            panel.hidden = panel.dataset.panel !== tab.dataset.view;
                        }});
                    }});
                }});

                updateFrames();
            }})();
        </script>
    </div>
    """


def _export_plotter_sequence_assets(
    plotters: _Sequence[_pv.Plotter],
    target_directory: _Path,
    prefix: str,
    *,
    export_3d_view: bool,
) -> tuple[list[_Path], list[_Path] | None]:
    """Export screenshots and optional interactive scenes for all frames."""
    target_directory.mkdir(parents=True, exist_ok=True)
    static_paths = []
    interactive_paths = [] if export_3d_view else None

    for index, plotter in enumerate(plotters):
        static_path = target_directory / f"{prefix}-{index}.png"
        plotter.screenshot(static_path)
        static_paths.append(static_path)

        if interactive_paths is not None:
            interactive_path = target_directory / f"{prefix}-{index}.html"
            plotter.export_html(interactive_path)
            interactive_paths.append(interactive_path)

    return static_paths, interactive_paths


def export_plotter_sequence(
    plotters: _Sequence[_pv.Plotter],
    file_name: _Path | str,
    *,
    labels: _Sequence[str] | None = None,
    export_3d_view: bool = True,
) -> _Path:
    """Export a sequence of PyVista plotters as a standalone HTML viewer.

    The viewer provides previous/next buttons and, by default, tabs for static
    screenshots and interactive 3D scenes. Assets are written to a sibling
    directory named ``<file stem>_files``.

    Args:
        plotters: One PyVista plotter for each frame.
        file_name: Path of the HTML viewer to create.
        labels: Optional label for each frame.
        export_3d_view: Whether to include interactive 3D scenes.

    Returns:
        The path of the created HTML viewer.
    """
    frame_labels = _validate_plotter_sequence(plotters, labels)
    viewer_path = _Path(file_name)
    if viewer_path.suffix == "":
        viewer_path = viewer_path.with_suffix(".html")
    elif viewer_path.suffix.lower() != ".html":
        raise ValueError("The plotter sequence file must have the '.html' extension.")

    viewer_path.parent.mkdir(parents=True, exist_ok=True)
    asset_directory = viewer_path.parent / f"{viewer_path.stem}_files"
    static_paths, interactive_paths = _export_plotter_sequence_assets(
        plotters,
        asset_directory,
        "frame",
        export_3d_view=export_3d_view,
    )
    static_sources = [
        path.relative_to(viewer_path.parent).as_posix() for path in static_paths
    ]
    interactive_sources = (
        [
            path.relative_to(viewer_path.parent).as_posix()
            for path in interactive_paths
        ]
        if interactive_paths is not None
        else None
    )
    viewer_path.write_text(
        _sequence_viewer_html(
            static_sources,
            frame_labels,
            interactive_sources=interactive_sources,
        ),
        encoding="utf-8",
    )
    return viewer_path


def show_plotter_sequence(
    plotters: _Sequence[_pv.Plotter],
    *,
    labels: _Sequence[str] | None = None,
    nbsphinx_export_3d_view: bool = True,
) -> None:
    """Show a sequence of PyVista scenes with previous/next controls.

    In Jupyter, the complete viewer is embedded in the notebook output. During
    an nbsphinx build, its assets are written to the documentation's static
    directory and the viewer is embedded as an iframe.

    Args:
        plotters: One PyVista plotter for each frame.
        labels: Optional label for each frame.
        nbsphinx_export_3d_view: Whether to include interactive 3D scenes. The
            static frame sequence is always included.
    """
    frame_labels = _validate_plotter_sequence(plotters, labels)
    sequence_uid = f"sequence-{_uuid.uuid4().hex}"

    if _is_nbsphinx():
        static_doc_path = _Path(_os.environ["PYVISTA_DOCS_STATIC"])
        viewer_path = static_doc_path / f"{sequence_uid}.html"
        export_plotter_sequence(
            plotters,
            viewer_path,
            labels=frame_labels,
            export_3d_view=nbsphinx_export_3d_view,
        )
        _display(
            _IFrame(
                src=f"../_static/pyvista/{viewer_path.name}",
                width="100%",
                height=650,
            )
        )
        return

    with _tempfile.TemporaryDirectory() as temporary_directory:
        temporary_path = _Path(temporary_directory)
        static_paths, interactive_paths = _export_plotter_sequence_assets(
            plotters,
            temporary_path,
            sequence_uid,
            export_3d_view=nbsphinx_export_3d_view,
        )

        static_sources = [
            "data:image/png;base64,"
            + _base64.b64encode(path.read_bytes()).decode("ascii")
            for path in static_paths
        ]
        interactive_sources = (
            [
                "data:text/html;base64,"
                + _base64.b64encode(path.read_bytes()).decode("ascii")
                for path in interactive_paths
            ]
            if interactive_paths is not None
            else None
        )

    _display(
        _HTML(
            _sequence_viewer_html(
                static_sources,
                frame_labels,
                interactive_sources=interactive_sources,
            )
        )
    )


def show_plotter(plotter: _pv.Plotter, *, nbsphinx_export_3d_view: bool = True) -> None:
    """Show a PyVista plotter.

    This function displays the plotter according to the current environment:
    - For local development, it will directly show the plotter.
    - For nbsphinx documentation, it will export a static image of the plotter
      and optionally an interactive 3D view that can be embedded in the
      documentation.

    Args:
        plotter: The PyVista plotter to show.
        nbsphinx_export_3d_view: For nbsphinx documentation, whether to
            export an interactive 3D view alongside the static image.
    """
    if not _is_nbsphinx():
        plotter.show()
    else:
        # Path where static documents are stored for the website.
        static_doc_path = _Path(_os.environ["PYVISTA_DOCS_STATIC"])

        # Get a unique identifier for the current plotter
        plotter_uid = str(_uuid.uuid4().hex)

        # Export a screenshot of the plotter
        plotter.screenshot(static_doc_path / f"{plotter_uid}.png")
        static_frame_html = f"""
            <img src="../_static/pyvista/{plotter_uid}.png"
                style="max-width:100%; border-radius:8px;">
            """

        if nbsphinx_export_3d_view:
            # Export a html representation of the plotter
            plotter.export_html(static_doc_path / f"{plotter_uid}.html")

            interactive_frame = _IFrame(
                src=f"../_static/pyvista/{plotter_uid}.html",
                width="100%",
                height=600,
            )
            tab = _widgets.Tab(
                children=[
                    _widgets.HTML(static_frame_html),
                    _widgets.HTML(interactive_frame._repr_html_()),
                ]
            )
            tab.set_title(0, "Static Scene")
            tab.set_title(1, "Interactive Scene")

            _display(tab)
        else:
            _display(_HTML(static_frame_html))
