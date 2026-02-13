"""
Report Generator for TRNspot
============================

Main class and functions for generating comprehensive HTML and PDF reports.
"""

import os
import glob
import base64
import json
from datetime import datetime
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field
from pathlib import Path

from .. import config


@dataclass
class ReportSection:
    """A section in the report."""

    title: str
    content: str
    section_id: str = ""
    figures: List[str] = field(default_factory=list)
    tables: List[Dict[str, Any]] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    subsections: List["ReportSection"] = field(default_factory=list)

    def __post_init__(self):
        if not self.section_id:
            self.section_id = self.title.lower().replace(" ", "-").replace("/", "-")


@dataclass
class ReportData:
    """Container for all report data."""

    title: str = "TRNspot Analysis Report"
    subtitle: str = ""
    timestamp: str = field(
        default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    sections: List[ReportSection] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


class ReportGenerator:
    """
    Generate comprehensive analysis reports in HTML and PDF formats.

    This class collects information from the pipeline run and generates
    professional reports with interactive navigation, plot galleries,
    and detailed analysis summaries.

    Parameters
    ----------
    output_dir : str
        Directory for output files
    title : str
        Report title

    Examples
    --------
    >>> generator = ReportGenerator(output_dir="results", title="My Analysis")
    >>> generator.add_section(create_data_summary_section(adata))
    >>> generator.add_section(create_settings_section())
    >>> generator.generate_html("report.html")
    >>> generator.generate_pdf("report.pdf")
    """

    def __init__(
        self,
        output_dir: str,
        title: str = "TRNspot Analysis Report",
        subtitle: str = "",
    ):
        self.output_dir = output_dir
        self.report_data = ReportData(
            title=title,
            subtitle=subtitle,
        )
        self._figures_dir = os.path.join(output_dir, "figures")

    def add_section(self, section: ReportSection) -> None:
        """Add a section to the report."""
        self.report_data.sections.append(section)

    def add_metadata(self, key: str, value: Any) -> None:
        """Add metadata to the report."""
        self.report_data.metadata[key] = value

    def _image_to_base64(self, image_path: str) -> str:
        """Convert image to base64 for embedding in HTML."""
        if not os.path.exists(image_path):
            return ""

        with open(image_path, "rb") as f:
            data = f.read()

        ext = Path(image_path).suffix.lower()
        mime_types = {
            ".png": "image/png",
            ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg",
            ".svg": "image/svg+xml",
            ".gif": "image/gif",
        }
        mime_type = mime_types.get(ext, "image/png")

        return f"data:{mime_type};base64,{base64.b64encode(data).decode()}"

    def _get_css(self) -> str:
        """Get CSS styles for the HTML report."""
        return """
        <style>
            :root {
                --primary-color: #2563eb;
                --primary-dark: #1d4ed8;
                --secondary-color: #64748b;
                --accent-color: #10b981;
                --warning-color: #f59e0b;
                --error-color: #ef4444;
                --background-color: #f8fafc;
                --card-background: #ffffff;
                --text-color: #1e293b;
                --text-light: #64748b;
                --border-color: #e2e8f0;
                --code-bg: #f1f5f9;
            }

            * {
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }

            html {
                scroll-behavior: smooth;
            }

            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
                line-height: 1.6;
                color: var(--text-color);
                background-color: var(--background-color);
            }

            .container {
                max-width: 1400px;
                margin: 0 auto;
                padding: 2rem;
            }

            /* Header */
            header {
                background: linear-gradient(135deg, var(--primary-color), var(--primary-dark));
                color: white;
                padding: 3rem 2rem;
                margin-bottom: 2rem;
                border-radius: 0.75rem;
                box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
            }

            header h1 {
                font-size: 2.5rem;
                margin-bottom: 0.5rem;
                font-weight: 700;
            }

            header .subtitle {
                font-size: 1.2rem;
                opacity: 0.9;
            }

            header .timestamp {
                font-size: 0.9rem;
                opacity: 0.7;
                margin-top: 1rem;
            }

            /* Navigation */
            .sidebar {
                position: fixed;
                left: 0;
                top: 0;
                width: 280px;
                height: 100vh;
                background: var(--card-background);
                border-right: 1px solid var(--border-color);
                padding: 1.5rem;
                overflow-y: auto;
                z-index: 1000;
            }

            .sidebar-header {
                font-size: 1.1rem;
                font-weight: 600;
                color: var(--primary-color);
                margin-bottom: 1rem;
                padding-bottom: 0.5rem;
                border-bottom: 2px solid var(--primary-color);
            }

            .sidebar ul {
                list-style: none;
            }

            .sidebar li {
                margin: 0.25rem 0;
            }

            .sidebar a {
                display: block;
                padding: 0.5rem 0.75rem;
                color: var(--text-color);
                text-decoration: none;
                border-radius: 0.375rem;
                transition: all 0.2s;
                font-size: 0.9rem;
            }

            .sidebar a:hover {
                background: var(--background-color);
                color: var(--primary-color);
            }

            .sidebar a.active {
                background: var(--primary-color);
                color: white;
            }

            .sidebar .subsection {
                padding-left: 1.5rem;
                font-size: 0.85rem;
            }

            .main-content {
                margin-left: 300px;
            }

            /* Sections */
            .section {
                background: var(--card-background);
                border-radius: 0.75rem;
                padding: 2rem;
                margin-bottom: 2rem;
                box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            }

            .section h2 {
                color: var(--primary-color);
                font-size: 1.5rem;
                margin-bottom: 1.5rem;
                padding-bottom: 0.75rem;
                border-bottom: 2px solid var(--border-color);
                display: flex;
                align-items: center;
                gap: 0.5rem;
            }

            .section h3 {
                font-size: 1.2rem;
                margin: 1.5rem 0 0.75rem;
                color: var(--text-color);
            }

            .section h4 {
                font-size: 1rem;
                margin: 1rem 0 0.5rem;
                color: var(--secondary-color);
            }

            /* Metrics Grid */
            .metrics-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 1rem;
                margin: 1.5rem 0;
            }

            .metric-card {
                background: linear-gradient(135deg, var(--background-color), #ffffff);
                padding: 1.25rem;
                border-radius: 0.5rem;
                text-align: center;
                border: 1px solid var(--border-color);
                transition: transform 0.2s;
            }

            .metric-card:hover {
                transform: translateY(-2px);
            }

            .metric-card .value {
                font-size: 2rem;
                font-weight: 700;
                color: var(--primary-color);
                line-height: 1.2;
            }

            .metric-card .label {
                font-size: 0.85rem;
                color: var(--secondary-color);
                margin-top: 0.25rem;
            }

            /* Status Badges */
            .status-badge {
                display: inline-flex;
                align-items: center;
                gap: 0.25rem;
                padding: 0.25rem 0.75rem;
                border-radius: 9999px;
                font-size: 0.8rem;
                font-weight: 500;
            }

            .status-completed {
                background: #dcfce7;
                color: #166534;
            }

            .status-skipped {
                background: #fef3c7;
                color: #92400e;
            }

            .status-error {
                background: #fee2e2;
                color: #991b1b;
            }

            /* Tables */
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 1rem 0;
                font-size: 0.9rem;
            }

            th, td {
                padding: 0.75rem 1rem;
                text-align: left;
                border-bottom: 1px solid var(--border-color);
            }

            th {
                background: var(--background-color);
                font-weight: 600;
                color: var(--text-color);
            }

            tr:hover {
                background-color: var(--background-color);
            }

            /* Code blocks */
            code {
                background: var(--code-bg);
                padding: 0.125rem 0.375rem;
                border-radius: 0.25rem;
                font-family: 'SF Mono', Monaco, 'Cascadia Code', monospace;
                font-size: 0.875rem;
            }

            pre {
                background: var(--code-bg);
                padding: 1rem;
                border-radius: 0.5rem;
                overflow-x: auto;
                font-size: 0.85rem;
            }

            pre code {
                background: none;
                padding: 0;
            }

            /* Figure Gallery - Compact Thumbnail Grid */
            .gallery {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(160px, 1fr));
                gap: 0.75rem;
                margin: 1rem 0;
                padding: 1rem;
                background: var(--background-color);
                border-radius: 0.5rem;
            }

            .gallery-item {
                background: var(--card-background);
                border: 1px solid var(--border-color);
                border-radius: 0.375rem;
                overflow: hidden;
                transition: all 0.2s ease;
                cursor: pointer;
                aspect-ratio: 1;
            }

            .gallery-item:hover {
                transform: scale(1.05);
                box-shadow: 0 8px 16px rgba(0, 0, 0, 0.15);
                z-index: 10;
                border-color: var(--primary-color);
            }

            .gallery-item img {
                width: 100%;
                height: 120px;
                object-fit: contain;
                background: #fafafa;
                padding: 0.25rem;
            }

            .gallery-item .caption {
                padding: 0.35rem 0.5rem;
                font-size: 0.7rem;
                color: var(--text-color);
                text-align: center;
                border-top: 1px solid var(--border-color);
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
                background: var(--card-background);
            }

            /* Gallery Category Header */
            .gallery-category {
                margin: 1.5rem 0 0.5rem 0;
            }

            .gallery-category-title {
                font-size: 0.9rem;
                font-weight: 600;
                color: var(--secondary-color);
                text-transform: uppercase;
                letter-spacing: 0.05em;
                padding-bottom: 0.5rem;
                border-bottom: 1px solid var(--border-color);
                margin-bottom: 0.5rem;
            }

            .gallery-category-count {
                font-size: 0.75rem;
                font-weight: 400;
                color: var(--text-light);
                margin-left: 0.5rem;
            }

            /* Larger gallery variant for main figures */
            .gallery-large {
                grid-template-columns: repeat(auto-fill, minmax(240px, 1fr));
                gap: 1rem;
            }

            .gallery-large .gallery-item img {
                height: 180px;
            }

            .gallery-large .gallery-item .caption {
                font-size: 0.8rem;
                padding: 0.5rem;
            }

            /* Subsection Tabs Navigation */
            .tabs-container {
                margin-top: 1.5rem;
            }

            .tabs-nav {
                display: flex;
                flex-wrap: wrap;
                gap: 0.25rem;
                border-bottom: 2px solid var(--border-color);
                margin-bottom: 0;
            }

            .tab-button {
                padding: 0.75rem 1.25rem;
                background: transparent;
                border: none;
                border-bottom: 3px solid transparent;
                cursor: pointer;
                font-size: 0.9rem;
                font-weight: 500;
                color: var(--secondary-color);
                transition: all 0.2s;
                margin-bottom: -2px;
            }

            .tab-button:hover {
                color: var(--primary-color);
                background: var(--background-color);
            }

            .tab-button.active {
                color: var(--primary-color);
                border-bottom-color: var(--primary-color);
                background: var(--card-background);
            }

            .tab-content {
                display: none;
                padding: 1.5rem;
                background: var(--background-color);
                border: 1px solid var(--border-color);
                border-top: none;
                border-radius: 0 0 0.5rem 0.5rem;
            }

            .tab-content.active {
                display: block;
            }

            .tab-content h3 {
                color: var(--primary-color);
                font-size: 1.1rem;
                margin: 0 0 1rem 0;
                padding-bottom: 0.5rem;
                border-bottom: 1px solid var(--border-color);
            }

            /* Legacy subsection styling (for non-tabbed content) */
            .subsection {
                background: var(--background-color);
                border: 1px solid var(--border-color);
                border-radius: 0.5rem;
                padding: 1.5rem;
                margin: 1rem 0;
            }

            .subsection h3 {
                color: var(--primary-color);
                font-size: 1.1rem;
                margin: 0 0 1rem 0;
                padding-bottom: 0.5rem;
                border-bottom: 1px solid var(--border-color);
            }

            /* Lightbox */
            .lightbox {
                display: none;
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background: rgba(0, 0, 0, 0.9);
                z-index: 2000;
                justify-content: center;
                align-items: center;
            }

            .lightbox.active {
                display: flex;
            }

            .lightbox img {
                max-width: 90%;
                max-height: 90%;
                object-fit: contain;
            }

            .lightbox-close {
                position: absolute;
                top: 1rem;
                right: 1rem;
                color: white;
                font-size: 2rem;
                cursor: pointer;
                background: none;
                border: none;
            }

            .lightbox-nav {
                position: absolute;
                top: 50%;
                transform: translateY(-50%);
                color: white;
                font-size: 3rem;
                cursor: pointer;
                background: none;
                border: none;
                padding: 1rem;
            }

            .lightbox-prev { left: 1rem; }
            .lightbox-next { right: 1rem; }

            /* Operations Log */
            .operations-log {
                max-height: 400px;
                overflow-y: auto;
                border: 1px solid var(--border-color);
                border-radius: 0.5rem;
            }

            .log-entry {
                display: flex;
                align-items: center;
                padding: 0.75rem 1rem;
                border-bottom: 1px solid var(--border-color);
                gap: 1rem;
            }

            .log-entry:last-child {
                border-bottom: none;
            }

            .log-entry .timestamp {
                font-family: monospace;
                font-size: 0.8rem;
                color: var(--secondary-color);
                white-space: nowrap;
            }

            .log-entry .step {
                font-weight: 500;
                flex: 1;
            }

            /* Settings Grid */
            .settings-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                gap: 1.5rem;
            }

            .settings-category {
                background: var(--background-color);
                padding: 1.25rem;
                border-radius: 0.5rem;
            }

            .settings-category h4 {
                color: var(--primary-color);
                margin-bottom: 0.75rem;
                font-size: 0.95rem;
            }

            .settings-category dl {
                display: grid;
                grid-template-columns: 1fr 1fr;
                gap: 0.5rem;
                font-size: 0.85rem;
            }

            .settings-category dt {
                color: var(--secondary-color);
            }

            .settings-category dd {
                font-weight: 500;
                text-align: right;
            }

            /* Footer */
            footer {
                text-align: center;
                padding: 2rem;
                color: var(--secondary-color);
                font-size: 0.875rem;
            }

            /* Responsive */
            @media (max-width: 1024px) {
                .sidebar {
                    display: none;
                }
                .main-content {
                    margin-left: 0;
                }
            }

            @media print {
                body { background: white; }
                .container { max-width: none; padding: 0; }
                .sidebar { display: none; }
                .main-content { margin-left: 0; }
                .section { box-shadow: none; page-break-inside: avoid; }
                header { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
                .gallery { grid-template-columns: repeat(2, 1fr); }
                .lightbox { display: none !important; }
            }
        </style>
        """

    def _get_javascript(self) -> str:
        """Get JavaScript for interactive features."""
        return """
        <script>
            // Lightbox functionality
            let currentGallery = [];
            let currentIndex = 0;

            function openLightbox(gallery, index) {
                currentGallery = gallery;
                currentIndex = index;
                updateLightboxImage();
                document.getElementById('lightbox').classList.add('active');
                document.body.style.overflow = 'hidden';
            }

            function closeLightbox() {
                document.getElementById('lightbox').classList.remove('active');
                document.body.style.overflow = '';
            }

            function navigateLightbox(direction) {
                currentIndex = (currentIndex + direction + currentGallery.length) % currentGallery.length;
                updateLightboxImage();
            }

            function updateLightboxImage() {
                document.getElementById('lightbox-img').src = currentGallery[currentIndex];
            }

            // Keyboard navigation
            document.addEventListener('keydown', function(e) {
                if (!document.getElementById('lightbox').classList.contains('active')) return;
                if (e.key === 'Escape') closeLightbox();
                if (e.key === 'ArrowLeft') navigateLightbox(-1);
                if (e.key === 'ArrowRight') navigateLightbox(1);
            });

            // Tab navigation functionality
            function switchTab(containerId, tabId) {
                const container = document.getElementById(containerId);
                if (!container) return;
                
                // Hide all tab contents and deactivate all buttons
                container.querySelectorAll('.tab-content').forEach(content => {
                    content.classList.remove('active');
                });
                container.querySelectorAll('.tab-button').forEach(btn => {
                    btn.classList.remove('active');
                });
                
                // Show selected tab content and activate button
                const tabContent = document.getElementById(tabId);
                const tabButton = container.querySelector(`[data-tab="${tabId}"]`);
                if (tabContent) tabContent.classList.add('active');
                if (tabButton) tabButton.classList.add('active');
            }

            // Active navigation highlighting
            const observer = new IntersectionObserver((entries) => {
                entries.forEach(entry => {
                    const id = entry.target.getAttribute('id');
                    const navLink = document.querySelector(`.sidebar a[href="#${id}"]`);
                    if (navLink) {
                        if (entry.isIntersecting) {
                            document.querySelectorAll('.sidebar a').forEach(a => a.classList.remove('active'));
                            navLink.classList.add('active');
                        }
                    }
                });
            }, { threshold: 0.1 });

            document.querySelectorAll('.section[id]').forEach(section => {
                observer.observe(section);
            });
        </script>
        """

    def _render_metrics(self, metrics: Dict[str, Any]) -> str:
        """Render metrics as HTML cards."""
        if not metrics:
            return ""

        cards = []
        for label, value in metrics.items():
            if isinstance(value, float):
                value_str = f"{value:,.2f}"
            elif isinstance(value, int):
                value_str = f"{value:,}"
            else:
                value_str = str(value)

            cards.append(
                f"""
            <div class="metric-card">
                <div class="value">{value_str}</div>
                <div class="label">{label}</div>
            </div>
            """
            )

        return f'<div class="metrics-grid">{"".join(cards)}</div>'

    def _render_figures(
        self, figures: List[str], embed: bool = True, gallery_id: str = "gallery"
    ) -> str:
        """Render figures as an interactive gallery.

        Images are stored once in a JavaScript array and referenced by index
        to avoid duplicating base64 data in each gallery item's onclick handler.
        """
        if not figures:
            return ""

        valid_figures = [f for f in figures if os.path.exists(f)]
        if not valid_figures:
            return "<p><em>No figures available</em></p>"

        # Pre-compute image sources once
        if embed:
            image_sources = [self._image_to_base64(f) for f in valid_figures]
        else:
            image_sources = valid_figures

        # Create a unique array name for this gallery to avoid conflicts
        array_name = f"galleryImages_{gallery_id.replace('-', '_')}"

        # Build JavaScript array with images (stored only once)
        images_json = json.dumps(image_sources)
        html_parts = [
            f"<script>var {array_name} = {images_json};</script>",
            f'<div class="gallery" id="{gallery_id}">',
        ]

        for idx, fig_path in enumerate(valid_figures):
            caption = Path(fig_path).stem.replace("_", " ").replace("-", " ").title()
            if len(caption) > 40:
                caption = caption[:37] + "..."

            # Reference the pre-stored array by index instead of embedding all images again
            gallery_js = f"openLightbox({array_name}, {idx})"

            html_parts.append(
                f"""
            <div class="gallery-item" onclick='{gallery_js}'>
                <img src="{image_sources[idx]}" alt="{caption}" loading="lazy">
                <div class="caption" title="{Path(fig_path).stem}">{caption}</div>
            </div>
            """
            )

        html_parts.append("</div>")
        return "".join(html_parts)

    def _render_tables(self, tables: List[Dict[str, Any]]) -> str:
        """Render tables as HTML."""
        if not tables:
            return ""

        html_parts = []
        for table_data in tables:
            title = table_data.get("title", "")
            headers = table_data.get("headers", [])
            rows = table_data.get("rows", [])

            if title:
                html_parts.append(f"<h4>{title}</h4>")

            header_html = "".join(f"<th>{h}</th>" for h in headers)
            rows_html = ""
            for row in rows:
                cells = "".join(f"<td>{cell}</td>" for cell in row)
                rows_html += f"<tr>{cells}</tr>"

            html_parts.append(
                f"""
            <table>
                <thead><tr>{header_html}</tr></thead>
                <tbody>{rows_html}</tbody>
            </table>
            """
            )

        return "".join(html_parts)

    def _render_section(self, section: ReportSection, idx: int) -> str:
        """Render a single section as HTML with tabbed navigation for subsections."""
        section_id = section.section_id or f"section-{idx}"

        subsections_html = ""
        if section.subsections:
            # Create tab navigation container
            tabs_container_id = f"{section_id}-tabs"

            # Build tab buttons
            tab_buttons = []
            for sub_idx, subsection in enumerate(section.subsections):
                is_active = "active" if sub_idx == 0 else ""
                tab_id = f"{section_id}-tab-{sub_idx}"
                tab_buttons.append(
                    f'<button class="tab-button {is_active}" data-tab="{tab_id}" '
                    f"onclick=\"switchTab('{tabs_container_id}', '{tab_id}')\">{subsection.title}</button>"
                )

            # Build tab contents
            tab_contents = []
            for sub_idx, subsection in enumerate(section.subsections):
                is_active = "active" if sub_idx == 0 else ""
                tab_id = f"{section_id}-tab-{sub_idx}"
                tab_contents.append(
                    f"""
                <div class="tab-content {is_active}" id="{tab_id}">
                    {self._render_metrics(subsection.metrics)}
                    <div class="content">{subsection.content}</div>
                    {self._render_figures(subsection.figures, gallery_id=f"{section_id}-gallery-{sub_idx}")}
                    {self._render_tables(subsection.tables)}
                </div>
                """
                )

            subsections_html = f"""
            <div class="tabs-container" id="{tabs_container_id}">
                <div class="tabs-nav">
                    {"".join(tab_buttons)}
                </div>
                {"".join(tab_contents)}
            </div>
            """

        return f"""
        <div class="section" id="{section_id}">
            <h2>{section.title}</h2>
            {self._render_metrics(section.metrics)}
            <div class="content">{section.content}</div>
            {self._render_figures(section.figures, gallery_id=f"{section_id}-gallery")}
            {self._render_tables(section.tables)}
            {subsections_html}
        </div>
        """

    def _render_sidebar(self) -> str:
        """Render navigation sidebar with main sections only (subsections are tabs)."""
        items = []
        for section in self.report_data.sections:
            items.append(
                f'<li><a href="#{section.section_id}">{section.title}</a></li>'
            )
            # Subsections are now rendered as tabs within their parent section,
            # so we don't add them to the sidebar navigation

        return f"""
        <nav class="sidebar">
            <div class="sidebar-header">Contents</div>
            <ul>{"".join(items)}</ul>
        </nav>
        """

    def _render_lightbox(self) -> str:
        """Render lightbox HTML for image viewer."""
        return """
        <div id="lightbox" class="lightbox" onclick="closeLightbox()">
            <button class="lightbox-close" onclick="closeLightbox()">&times;</button>
            <button class="lightbox-nav lightbox-prev" onclick="event.stopPropagation(); navigateLightbox(-1)">&#8249;</button>
            <img id="lightbox-img" src="" alt="Full size image" onclick="event.stopPropagation()">
            <button class="lightbox-nav lightbox-next" onclick="event.stopPropagation(); navigateLightbox(1)">&#8250;</button>
        </div>
        """

    def generate_html(
        self,
        output_path: Optional[str] = None,
        embed_images: bool = True,
    ) -> str:
        """
        Generate HTML report.

        Parameters
        ----------
        output_path : str, optional
            Path to save the HTML file
        embed_images : bool
            Whether to embed images as base64

        Returns
        -------
        str
            HTML content
        """
        sections_html = "\n".join(
            self._render_section(section, idx)
            for idx, section in enumerate(self.report_data.sections)
        )

        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{self.report_data.title}</title>
            {self._get_css()}
        </head>
        <body>
            {self._render_sidebar()}
            
            <div class="main-content">
                <div class="container">
                    <header>
                        <h1>{self.report_data.title}</h1>
                        <div class="subtitle">{self.report_data.subtitle}</div>
                        <div class="timestamp">Generated: {self.report_data.timestamp}</div>
                    </header>

                    {sections_html}

                    <footer>
                        Generated by TRNspot v{getattr(config, '__version__', '0.1.0')} | 
                        Transcriptional Regulatory Network Analysis
                    </footer>
                </div>
            </div>

            {self._render_lightbox()}
            {self._get_javascript()}
        </body>
        </html>
        """

        if output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html)
            print(f"HTML report saved to: {output_path}")

        return html

    def generate_pdf(
        self,
        output_path: Optional[str] = None,
    ) -> Optional[str]:
        """
        Generate PDF report from HTML.

        Requires weasyprint or pdfkit to be installed.

        Parameters
        ----------
        output_path : str, optional
            Path to save the PDF file

        Returns
        -------
        str or None
            Path to PDF file, or None if generation failed
        """
        if output_path is None:
            output_path = os.path.join(self.output_dir, "report.pdf")

        html_content = self.generate_html(embed_images=True)

        # Try weasyprint first
        try:
            from weasyprint import HTML

            HTML(string=html_content).write_pdf(output_path)
            print(f"✓ PDF report saved to: {output_path}")
            return output_path
        except ImportError:
            pass
        except Exception as e:
            print(f"Warning: weasyprint PDF generation failed: {e}")

        # Try pdfkit as fallback
        try:
            import pdfkit

            pdfkit.from_string(html_content, output_path)
            print(f"✓ PDF report saved to: {output_path}")
            return output_path
        except ImportError:
            pass
        except Exception as e:
            print(f"Warning: pdfkit PDF generation failed: {e}")

        print("⚠ PDF generation requires 'weasyprint' or 'pdfkit'. Install with:")
        print("  pip install weasyprint")
        print("  or")
        print("  pip install pdfkit")
        return None


def generate_report(
    output_dir: str,
    title: str = "TRNspot Analysis Report",
    subtitle: str = "",
    adata=None,
    celloracle_result=None,
    hotspot_result=None,
    merged_scores=None,
    log_file: Optional[str] = None,
    formats: List[str] = ["html", "pdf"],
    embed_images: bool = True,
) -> Dict[str, str]:
    """
    Generate comprehensive analysis report in specified formats.

    Parameters
    ----------
    output_dir : str
        Output directory
    title : str
        Report title
    subtitle : str
        Report subtitle
    adata : AnnData, optional
        Processed AnnData object
    celloracle_result : tuple, optional
        CellOracle results (oracle, links)
    hotspot_result : Hotspot, optional
        Hotspot results
    merged_scores : pd.DataFrame, optional
        Merged GRN scores from cross-stratification analysis
    log_file : str, optional
        Path to pipeline log file
    formats : list
        Output formats ('html', 'pdf')
    embed_images : bool
        If True, embed images as base64 (larger file, self-contained).
        If False, use relative paths (smaller file, requires images in place).

    Returns
    -------
    dict
        Paths to generated report files
    """
    from .sections import (
        create_data_summary_section,
        create_settings_section,
        create_qc_section,
        create_preprocessing_section,
        create_clustering_section,
        create_celloracle_section,
        create_hotspot_section,
        create_operations_log_section,
        create_plot_gallery_section,
        create_grn_deep_analysis_section,
    )

    generator = ReportGenerator(output_dir=output_dir, title=title, subtitle=subtitle)

    # Add sections based on available data
    if adata is not None:
        generator.add_section(create_data_summary_section(adata, output_dir))

    generator.add_section(create_settings_section())

    if adata is not None:
        generator.add_section(create_qc_section(adata, output_dir))
        generator.add_section(create_preprocessing_section(adata, output_dir))
        generator.add_section(create_clustering_section(adata, output_dir))

    if celloracle_result is not None:
        generator.add_section(create_celloracle_section(celloracle_result, output_dir))

    if hotspot_result is not None:
        generator.add_section(create_hotspot_section(hotspot_result, output_dir))

    # Add GRN deep analysis section if merged scores or heatmap figures exist
    merged_csv = os.path.join(output_dir, "celloracle", "total_merged_scores.csv")
    grn_heatmaps_exist = (
        any(
            "grn_heatmap_" in f
            for f in glob.glob(
                os.path.join(output_dir, "figures", "grn", "**", "*.png"),
                recursive=True,
            )
        )
        if os.path.exists(os.path.join(output_dir, "figures", "grn"))
        else False
    )

    if merged_scores is not None or os.path.exists(merged_csv) or grn_heatmaps_exist:
        generator.add_section(
            create_grn_deep_analysis_section(output_dir, merged_scores)
        )

    if log_file and os.path.exists(log_file):
        generator.add_section(create_operations_log_section(log_file))

    generator.add_section(create_plot_gallery_section(output_dir))

    # Generate reports
    outputs = {}

    if "html" in formats:
        html_path = os.path.join(output_dir, "report.html")
        generator.generate_html(html_path, embed_images=embed_images)
        outputs["html"] = html_path

    if "pdf" in formats:
        pdf_path = os.path.join(output_dir, "report.pdf")
        result = generator.generate_pdf(pdf_path)
        if result:
            outputs["pdf"] = pdf_path

    return outputs


def generate_html_report(output_dir: str, embed_images: bool = True, **kwargs) -> str:
    """Convenience function to generate HTML report only."""
    outputs = generate_report(
        output_dir, formats=["html"], embed_images=embed_images, **kwargs
    )
    return outputs.get("html", "")


def generate_pdf_report(output_dir: str, **kwargs) -> Optional[str]:
    """Convenience function to generate PDF report only."""
    outputs = generate_report(output_dir, formats=["pdf"], **kwargs)
    return outputs.get("pdf")
