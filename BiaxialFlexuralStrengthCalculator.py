"""
Biaxial Flexural Strength Calculator
=====================================
Calculates the biaxial flexural strength of ceramic disc specimens
using the closed-form stress solutions for ring-on-ring (ASTM C1499),
piston-on-three-balls (ISO 6872) and ball-on-three-balls (B3B) test
configurations. Includes two-parameter Weibull statistical analysis
by maximum likelihood estimation for batch data.

Author : Dr Thanos Goulas
Contact: thanosgoulas@outlook.com
Version: 1.0
Date   : 2026

Feedback, bug reports and suggestions are welcome. Please open an
issue on the project repository or get in touch by email.

SPDX-License-Identifier: MIT

MIT License

Copyright (c) 2026 Dr Thanos Goulas

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Disclaimer
----------
This software is provided as an engineering aid. The accuracy of
calculated values depends on the validity of the input parameters and
the underlying assumptions of the chosen test standard. The user is
responsible for verifying results against the relevant standard
(ASTM C1499, ISO 6872, etc.) and for confirming the appropriateness
of the test method and statistical treatment for their application.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import math
import csv
import json
from pathlib import Path

try:
    import numpy as np
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False


# ============================================================
# Windows 11 inspired theme
# ============================================================
THEME = {
    "bg":            "#F3F3F3",   # Mica-like background
    "card":          "#FBFBFB",   # Card / panel
    "card_border":   "#E5E5E5",
    "accent":        "#0067C0",   # Win11 system accent (blue)
    "accent_hover":  "#1976D2",
    "accent_press":  "#005A9E",
    "text":          "#1B1B1B",
    "text_muted":    "#5C5C5C",
    "text_subtle":   "#8A8A8A",
    "input_bg":      "#FFFFFF",
    "input_border":  "#D1D1D1",
    "input_focus":   "#0067C0",
    "success":       "#107C10",
    "warning":       "#C29200",
    "error":         "#C42B1C",
    "divider":       "#EBEBEB",
    "tooltip_bg":    "#2B2B2B",
    "tooltip_fg":    "#FFFFFF",
}

FONT_FAMILY = "Segoe UI Variable"  # Falls back to Segoe UI on older systems
FONT_FALLBACK = "Segoe UI"


# ============================================================
# Tooltip widget
# ============================================================
class Tooltip:
    """Tooltip widget. Displays a dark popup on hover after a short delay."""
    def __init__(self, widget, text, delay=450, wraplength=320):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wraplength = wraplength
        self.tip = None
        self.after_id = None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._hide)
        widget.bind("<ButtonPress>", self._hide)

    def _schedule(self, _=None):
        self._cancel()
        self.after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self.after_id:
            self.widget.after_cancel(self.after_id)
            self.after_id = None

    def _show(self):
        if self.tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 16
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{x}+{y}")
        self.tip.attributes("-topmost", True)
        frame = tk.Frame(self.tip, bg=THEME["tooltip_bg"], bd=0,
                         highlightthickness=1, highlightbackground="#444")
        frame.pack()
        label = tk.Label(frame, text=self.text, justify="left",
                         bg=THEME["tooltip_bg"], fg=THEME["tooltip_fg"],
                         font=(FONT_FALLBACK, 9), padx=10, pady=6,
                         wraplength=self.wraplength)
        label.pack()

    def _hide(self, _=None):
        self._cancel()
        if self.tip:
            self.tip.destroy()
            self.tip = None


# ============================================================
# Calculation engines
# ============================================================
class StrengthCalculator:
    """Closed-form solutions for biaxial flexural strength tests."""

    @staticmethod
    def ring_on_ring(F, h, nu, D_S, D_L, D):
        """
        ASTM C1499 ring-on-ring or piston-on-ring.
        sigma = (3F / 2 pi h^2) * [(1-nu)(D_S^2 - D_L^2)/(2 D^2) + (1+nu) ln(D_S/D_L)]
        Inputs in SI: F [N], lengths [m]. Returns Pa.
        """
        if h <= 0 or D <= 0 or D_L <= 0 or D_S <= 0:
            raise ValueError("All dimensions must be positive.")
        if D_L >= D_S:
            raise ValueError("Load ring diameter D_L must be smaller than support ring diameter D_S.")
        if D_S > D:
            raise ValueError("Support ring diameter D_S cannot exceed specimen diameter D.")
        term1 = (1 - nu) * (D_S**2 - D_L**2) / (2 * D**2)
        term2 = (1 + nu) * math.log(D_S / D_L)
        sigma = (3 * F) / (2 * math.pi * h**2) * (term1 + term2)
        return sigma

    @staticmethod
    def piston_on_three_balls(F, h, nu, r_support, r_load, r_specimen):
        """
        ISO 6872 piston-on-three-balls (Wachtman / dental ceramic standard).
        sigma = -0.2387 * F * (X - Y) / h^2
        with
          X = (1 + nu) ln((r2/r3)^2) + (1 - nu) (r2^2 - r3^2) / (2 r1^2)
          Y = (1 + nu) [1 + ln((r1/r3)^2)] + (1 - nu) (r1^2 - r3^2) / r3^2
        where:
          r1 = radius of support circle (centres of three balls)
          r2 = radius of loaded area at tensile surface; ISO 6872 recommends
               r2 = h/3 when the piston contact radius is unknown
          r3 = radius of specimen
        Inputs in SI: F [N], lengths [m]. Returns Pa.
        """
        if h <= 0 or r_support <= 0 or r_load <= 0 or r_specimen <= 0:
            raise ValueError("All dimensions must be positive.")
        if r_support >= r_specimen:
            raise ValueError("Support radius r1 must be smaller than specimen radius r3.")
        r1, r2, r3 = r_support, r_load, r_specimen
        X = (1 + nu) * math.log((r2 / r3)**2) + (1 - nu) * (r2**2 - r3**2) / (2 * r1**2)
        Y = (1 + nu) * (1 + math.log((r1 / r3)**2)) + (1 - nu) * (r1**2 - r3**2) / r3**2
        sigma = -0.2387 * F * (X - Y) / h**2
        return sigma

    @staticmethod
    def ball_on_three_balls(F, h, f_factor):
        """
        B3B test (Börger, Supancic, Danzer).
        sigma = f * F / h^2
        f is a dimensionless factor obtained from FEA based on
        specimen radius, support ball radius and Poisson's ratio.
        Inputs: F [N], h [m], f [dimensionless].
        Returns Pa.
        """
        if h <= 0:
            raise ValueError("Thickness must be positive.")
        if f_factor <= 0:
            raise ValueError("f-factor must be positive.")
        return f_factor * F / h**2


# ============================================================
# Weibull statistics
# ============================================================
class WeibullAnalysis:
    """Two-parameter Weibull analysis using maximum likelihood estimation."""

    @staticmethod
    def fit(strengths):
        """
        Maximum likelihood estimation of Weibull modulus m and characteristic
        strength sigma_0 for a 2-parameter distribution.
        Returns (m, sigma_0, n).
        """
        if not HAS_PLOTTING:
            raise RuntimeError("NumPy is required for Weibull analysis.")
        x = np.asarray([s for s in strengths if s > 0], dtype=float)
        n = len(x)
        if n < 2:
            raise ValueError("Need at least 2 valid strength values.")

        # Newton-Raphson on the MLE equation for m
        lnx = np.log(x)
        m = 1.0  # initial guess
        for _ in range(200):
            xm = x ** m
            S = xm.sum()
            S1 = (xm * lnx).sum()
            S2 = (xm * lnx * lnx).sum()
            f = S1 / S - lnx.mean() - 1.0 / m
            fprime = (S2 * S - S1 * S1) / (S * S) + 1.0 / (m * m)
            step = f / fprime
            m_new = m - step
            if m_new <= 0:
                m_new = m / 2
            if abs(m_new - m) < 1e-9:
                m = m_new
                break
            m = m_new

        sigma_0 = (np.mean(x ** m)) ** (1.0 / m)
        return float(m), float(sigma_0), n

    @staticmethod
    def probability_of_failure(strengths):
        """Median rank probability of failure: P_i = (i - 0.5) / n."""
        x = np.asarray(sorted(strengths), dtype=float)
        n = len(x)
        ranks = np.arange(1, n + 1)
        P = (ranks - 0.5) / n
        return x, P


# ============================================================
# Custom widgets in Win11 style
# ============================================================
class WinButton(tk.Canvas):
    """Canvas-drawn button with hover and press states."""
    def __init__(self, parent, text, command=None, style="accent",
                 width=140, height=34, **kwargs):
        bg = kwargs.pop("bg", THEME["card"])
        super().__init__(parent, width=width, height=height,
                         bg=bg, highlightthickness=0, bd=0, **kwargs)
        self.text = text
        self.command = command
        self.style = style
        self.w = width
        self.h = height
        self.state = "normal"
        self.bind("<Enter>", lambda e: self._set_state("hover"))
        self.bind("<Leave>", lambda e: self._set_state("normal"))
        self.bind("<ButtonPress-1>", lambda e: self._set_state("press"))
        self.bind("<ButtonRelease-1>", self._on_release)
        self._draw()

    def _colors(self):
        if self.style == "accent":
            base = {
                "normal": (THEME["accent"], "#FFFFFF", THEME["accent"]),
                "hover":  (THEME["accent_hover"], "#FFFFFF", THEME["accent_hover"]),
                "press":  (THEME["accent_press"], "#E0E0E0", THEME["accent_press"]),
            }
        else:  # subtle / secondary
            base = {
                "normal": ("#FFFFFF", THEME["text"], THEME["input_border"]),
                "hover":  ("#F5F5F5", THEME["text"], "#B5B5B5"),
                "press":  ("#EAEAEA", THEME["text"], "#A0A0A0"),
            }
        return base[self.state]

    def _draw(self):
        self.delete("all")
        fill, fg, border = self._colors()
        r = 6
        self._round_rect(1, 1, self.w - 1, self.h - 1, r,
                         fill=fill, outline=border)
        self.create_text(self.w / 2, self.h / 2, text=self.text,
                         fill=fg, font=(FONT_FALLBACK, 10))

    def _round_rect(self, x1, y1, x2, y2, r, **kwargs):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1, x2, y1 + r,
            x2, y2 - r, x2, y2, x2 - r, y2, x1 + r, y2,
            x1, y2, x1, y2 - r, x1, y1 + r, x1, y1,
        ]
        return self.create_polygon(pts, smooth=True, **kwargs)

    def _set_state(self, s):
        self.state = s
        self._draw()

    def _on_release(self, event):
        self._set_state("hover")
        if 0 <= event.x <= self.w and 0 <= event.y <= self.h and self.command:
            self.command()


# ============================================================
# Main application
# ============================================================
class BiaxialApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Biaxial Flexural Strength Calculator")
        # Hide the default Tk feather icon
        try:
            blank = tk.PhotoImage(width=1, height=1)
            self.iconphoto(True, blank)
            self._blank_icon = blank  # keep reference so it isn't garbage-collected
        except Exception:
            pass
        self.geometry("1180x780")
        self.minsize(1000, 700)
        self.configure(bg=THEME["bg"])

        # State
        self.test_method = tk.StringVar(value="Piston-on-three-balls (ISO 6872)")
        self.unit_force = tk.StringVar(value="N")
        self.unit_length = tk.StringVar(value="mm")
        self.specimens = []  # list of dicts: {force, sigma_MPa, valid}

        self._setup_styles()
        self._build_ui()

    # --------------------------------------------------------
    def _setup_styles(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure("TCombobox",
                        fieldbackground=THEME["input_bg"],
                        background=THEME["input_bg"],
                        foreground=THEME["text"],
                        bordercolor=THEME["input_border"],
                        lightcolor=THEME["input_border"],
                        darkcolor=THEME["input_border"],
                        arrowcolor=THEME["text_muted"],
                        padding=6)
        style.map("TCombobox",
                  fieldbackground=[("readonly", THEME["input_bg"])],
                  foreground=[("readonly", THEME["text"])])

        style.configure("Win.Treeview",
                        background=THEME["input_bg"],
                        fieldbackground=THEME["input_bg"],
                        foreground=THEME["text"],
                        rowheight=26,
                        font=(FONT_FALLBACK, 9),
                        bordercolor=THEME["card_border"])
        style.configure("Win.Treeview.Heading",
                        background=THEME["card"],
                        foreground=THEME["text_muted"],
                        font=(FONT_FALLBACK, 9, "bold"),
                        relief="flat",
                        padding=6)
        style.map("Win.Treeview",
                  background=[("selected", "#CCE4F7")],
                  foreground=[("selected", THEME["text"])])
        style.map("Win.Treeview.Heading",
                  background=[("active", THEME["divider"])])

        style.configure("Win.TNotebook", background=THEME["bg"], borderwidth=0)
        style.configure("Win.TNotebook.Tab",
                        background=THEME["bg"],
                        foreground=THEME["text_muted"],
                        padding=(16, 8),
                        font=(FONT_FALLBACK, 10))
        style.map("Win.TNotebook.Tab",
                  background=[("selected", THEME["card"])],
                  foreground=[("selected", THEME["accent"])])

    # --------------------------------------------------------
    def _build_ui(self):
        # Title bar strip
        header = tk.Frame(self, bg=THEME["bg"], height=64)
        header.pack(fill="x", padx=20, pady=(16, 0))
        tk.Label(header,
                 text="Biaxial Flexural Strength Calculator",
                 bg=THEME["bg"], fg=THEME["text"],
                 font=(FONT_FALLBACK, 18, "bold")).pack(anchor="w")
        tk.Label(header,
                 text="ASTM C1499  ·  ISO 6872  ·  Ball-on-three-balls (B3B)  ·  Weibull statistics",
                 bg=THEME["bg"], fg=THEME["text_subtle"],
                 font=(FONT_FALLBACK, 10)).pack(anchor="w", pady=(2, 0))

        # Notebook
        notebook = ttk.Notebook(self, style="Win.TNotebook")
        notebook.pack(fill="both", expand=True, padx=20, pady=16)

        self.tab_calc = tk.Frame(notebook, bg=THEME["bg"])
        self.tab_weibull = tk.Frame(notebook, bg=THEME["bg"])
        self.tab_formulas = tk.Frame(notebook, bg=THEME["bg"])
        self.tab_about = tk.Frame(notebook, bg=THEME["bg"])
        notebook.add(self.tab_calc, text="  Single-specimen calculator  ")
        notebook.add(self.tab_weibull, text="  Batch & Weibull analysis  ")
        notebook.add(self.tab_formulas, text="  Formulas & references  ")
        notebook.add(self.tab_about, text="  About  ")

        self._build_calc_tab()
        self._build_weibull_tab()
        self._build_formulas_tab()
        self._build_about_tab()

    # --------------------------------------------------------
    # Helper: card frame
    def _card(self, parent, title=None):
        outer = tk.Frame(parent, bg=THEME["card_border"])
        inner = tk.Frame(outer, bg=THEME["card"])
        inner.pack(fill="both", expand=True, padx=1, pady=1)
        if title:
            tk.Label(inner, text=title, bg=THEME["card"], fg=THEME["text"],
                     font=(FONT_FALLBACK, 11, "bold")
                     ).pack(anchor="w", padx=18, pady=(14, 2))
        return outer, inner

    def _labeled_entry(self, parent, label, default="", tooltip=None,
                       width=14, info=True):
        row = tk.Frame(parent, bg=THEME["card"])
        row.pack(fill="x", padx=18, pady=4)
        lbl_frame = tk.Frame(row, bg=THEME["card"])
        lbl_frame.pack(side="left", fill="x", expand=True)
        lbl = tk.Label(lbl_frame, text=label,
                       bg=THEME["card"], fg=THEME["text"],
                       font=(FONT_FALLBACK, 10), anchor="w")
        lbl.pack(side="left")
        if tooltip and info:
            info_lbl = tk.Label(lbl_frame, text=" ⓘ",
                                bg=THEME["card"], fg=THEME["text_subtle"],
                                font=(FONT_FALLBACK, 10), cursor="hand2")
            info_lbl.pack(side="left")
            Tooltip(info_lbl, tooltip)
            Tooltip(lbl, tooltip)
        var = tk.StringVar(value=str(default))
        entry = tk.Entry(row, textvariable=var, width=width,
                         bg=THEME["input_bg"], fg=THEME["text"],
                         relief="solid", bd=1,
                         highlightthickness=1,
                         highlightbackground=THEME["input_border"],
                         highlightcolor=THEME["input_focus"],
                         insertbackground=THEME["text"],
                         font=(FONT_FALLBACK, 10))
        entry.pack(side="right", ipady=4)
        if tooltip:
            Tooltip(entry, tooltip)
        return var, entry

    # --------------------------------------------------------
    def _build_calc_tab(self):
        # Two columns: left = inputs, right = output + diagram
        container = tk.Frame(self.tab_calc, bg=THEME["bg"])
        container.pack(fill="both", expand=True)
        container.grid_columnconfigure(0, weight=1, uniform="cols")
        container.grid_columnconfigure(1, weight=1, uniform="cols")
        container.grid_rowconfigure(0, weight=1)

        left = tk.Frame(container, bg=THEME["bg"])
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        right = tk.Frame(container, bg=THEME["bg"])
        right.grid(row=0, column=1, sticky="nsew", padx=(10, 0))

        # --- Test method card ---
        method_card_outer, method_card = self._card(left, "Test configuration")
        method_card_outer.pack(fill="x", pady=(0, 12))

        row = tk.Frame(method_card, bg=THEME["card"])
        row.pack(fill="x", padx=18, pady=(0, 14))
        tk.Label(row, text="Test method",
                 bg=THEME["card"], fg=THEME["text"],
                 font=(FONT_FALLBACK, 10)).pack(side="left")
        info = tk.Label(row, text=" ⓘ",
                        bg=THEME["card"], fg=THEME["text_subtle"],
                        font=(FONT_FALLBACK, 10), cursor="hand2")
        info.pack(side="left")
        Tooltip(info,
                "Ring-on-Ring (ASTM C1499): equibiaxial stress state in the "
                "central load-ring region. Requires flat, parallel disc faces.\n\n"
                "Piston-on-three-balls (ISO 6872): adopted by ISO 6872 for "
                "dental ceramics. Less sensitive to specimen non-flatness than "
                "ring-on-ring.\n\n"
                "Ball-on-three-balls (B3B): suitable for small or as-sintered "
                "discs. Insensitive to thickness variation and surface waviness.")
        combo = ttk.Combobox(row, textvariable=self.test_method,
                             state="readonly", width=38,
                             values=[
                                 "Ring-on-Ring (ASTM C1499)",
                                 "Piston-on-three-balls (ISO 6872)",
                                 "Ball-on-three-balls (B3B)",
                             ])
        combo.pack(side="right")
        combo.bind("<<ComboboxSelected>>", lambda e: self._refresh_inputs())

        # --- Inputs card ---
        self.inputs_card_outer, self.inputs_card = self._card(left, "Specimen & loading parameters")
        self.inputs_card_outer.pack(fill="x", pady=(0, 12))

        # Units row
        units_row = tk.Frame(self.inputs_card, bg=THEME["card"])
        units_row.pack(fill="x", padx=18, pady=(0, 8))
        tk.Label(units_row, text="Units:", bg=THEME["card"],
                 fg=THEME["text_muted"],
                 font=(FONT_FALLBACK, 9)).pack(side="left")
        tk.Label(units_row, text="Force", bg=THEME["card"],
                 fg=THEME["text_muted"],
                 font=(FONT_FALLBACK, 9)).pack(side="left", padx=(10, 4))
        ttk.Combobox(units_row, textvariable=self.unit_force,
                     values=["N", "kN"], width=5,
                     state="readonly").pack(side="left")
        tk.Label(units_row, text="Length", bg=THEME["card"],
                 fg=THEME["text_muted"],
                 font=(FONT_FALLBACK, 9)).pack(side="left", padx=(10, 4))
        ttk.Combobox(units_row, textvariable=self.unit_length,
                     values=["mm", "m"], width=5,
                     state="readonly").pack(side="left")
        tk.Label(units_row, text="(stress output: MPa)",
                 bg=THEME["card"], fg=THEME["text_subtle"],
                 font=(FONT_FALLBACK, 9, "italic")).pack(side="left", padx=10)

        self.input_vars = {}  # populated by _refresh_inputs
        self.input_widgets_frame = tk.Frame(self.inputs_card, bg=THEME["card"])
        self.input_widgets_frame.pack(fill="x", pady=(4, 14))

        # --- Action buttons ---
        btn_row = tk.Frame(left, bg=THEME["bg"])
        btn_row.pack(fill="x", pady=(0, 12))
        WinButton(btn_row, "Calculate", command=self._calculate,
                  style="accent", width=130, bg=THEME["bg"]).pack(side="left")
        WinButton(btn_row, "Add to batch", command=self._add_to_batch,
                  style="subtle", width=120, bg=THEME["bg"]).pack(side="left", padx=6)
        WinButton(btn_row, "Save setup", command=self._save_setup,
                  style="subtle", width=110, bg=THEME["bg"]).pack(side="left", padx=(0, 6))
        WinButton(btn_row, "Load setup", command=self._load_setup,
                  style="subtle", width=110, bg=THEME["bg"]).pack(side="left")
        WinButton(btn_row, "Reset", command=lambda: self._reset_inputs(),
                  style="subtle", width=80, bg=THEME["bg"]).pack(side="left", padx=(6, 0))

        # --- Result card ---
        result_card_outer, result_card = self._card(right, "Result")
        result_card_outer.pack(fill="x", pady=(0, 12))
        self.result_value = tk.Label(result_card, text="—",
                                     bg=THEME["card"], fg=THEME["accent"],
                                     font=(FONT_FALLBACK, 32, "bold"))
        self.result_value.pack(anchor="w", padx=18, pady=(4, 0))
        self.result_unit = tk.Label(result_card, text="MPa",
                                    bg=THEME["card"], fg=THEME["text_muted"],
                                    font=(FONT_FALLBACK, 11))
        self.result_unit.pack(anchor="w", padx=18, pady=(0, 14))

        self.result_detail = tk.Label(result_card, text="",
                                      bg=THEME["card"], fg=THEME["text_muted"],
                                      font=(FONT_FALLBACK, 9), justify="left",
                                      wraplength=440, anchor="w")
        self.result_detail.pack(anchor="w", padx=18, pady=(0, 14), fill="x")

        # --- Schematic card ---
        schem_outer, schem_inner = self._card(right, "Test schematic")
        schem_outer.pack(fill="both", expand=True)
        self.schematic = tk.Canvas(schem_inner, bg=THEME["card"],
                                   highlightthickness=0,
                                   width=460, height=280)
        self.schematic.pack(fill="both", expand=True, padx=18, pady=(4, 14))
        # Redraw whenever the canvas is resized so the diagram scales to fit
        self.schematic.bind("<Configure>", lambda e: self._draw_schematic())

        self._refresh_inputs()

    # --------------------------------------------------------
    def _refresh_inputs(self, preserve=False):
        # Optionally preserve current values across rebuilds
        old_values = {}
        if preserve:
            for k, v in self.input_vars.items():
                try:
                    old_values[k] = v.get()
                except Exception:
                    pass
        for w in self.input_widgets_frame.winfo_children():
            w.destroy()
        self.input_vars = {}

        method = self.test_method.get()
        if method.startswith("Ring-on-Ring"):
            fields = [
                ("F",  "Failure load, F",       "1500",
                 "Maximum recorded force at failure.", True),
                ("h",  "Specimen thickness, h", "1.5",
                 "Disc thickness in the loaded region.", True),
                ("nu", "Poisson's ratio, ν",   "0.27",
                 "Poisson's ratio of the ceramic. Typical values: 3YSZ/4YSZ 0.28, "
                 "alumina 0.21, lithium disilicate 0.24, ATZ 0.29, ZTA 0.25.", False),
                ("DL", "Load ring diameter, D_L", "10.0",
                 "Diameter of the inner (loading) ring.", True),
                ("DS", "Support ring diameter, D_S", "20.0",
                 "Diameter of the outer (support) ring.", True),
                ("D",  "Specimen diameter, D",  "25.0",
                 "Disc diameter. Must be ≥ support ring diameter.", True),
            ]
        elif method.startswith("Piston"):
            fields = [
                ("F",  "Failure load, F",       "200",
                 "Maximum recorded force at failure.", True),
                ("h",  "Specimen thickness, h", "1.2",
                 "Disc thickness.", True),
                ("nu", "Poisson's ratio, ν",   "0.28",
                 "Poisson's ratio of the material. Typical values: "
                 "3YSZ/4YSZ 0.28, alumina 0.21, lithium disilicate 0.24, "
                 "ATZ 0.29, ZTA 0.25.", False),
                ("r1", "Support radius, r1 (preset)",   "6.2",
                 "Radius of the circle through the centres of the three support balls. "
                 "Default value 6.2 mm. Modify if a different fixture is used.", True),
                ("r2", "Piston radius, r2 (preset)", "0.8",
                 "Contact radius of the piston tip on the disc surface. Default value "
                 "0.8 mm (piston diameter 1.6 mm). Modify if a different fixture is "
                 "used.", True),
                ("D",  "Specimen diameter, D",  "13.0",
                 "Disc diameter. Converted to specimen radius (r3 = D/2) internally.", True),
            ]
        else:  # B3B
            fields = [
                ("F",  "Failure load, F",       "300",
                 "Maximum recorded force at failure.", True),
                ("h",  "Specimen thickness, h", "1.0",
                 "Disc thickness.", True),
                ("f",  "Geometry factor, f",    "1.55",
                 "Dimensionless factor from FEA (Börger, Supancic, Danzer). "
                 "Depends on specimen radius / support ball radius / Poisson's "
                 "ratio. Typical range 1.2–2.0 for common geometries.", False),
            ]

        # Two columns
        col_frame = tk.Frame(self.input_widgets_frame, bg=THEME["card"])
        col_frame.pack(fill="x")
        left_col = tk.Frame(col_frame, bg=THEME["card"])
        right_col = tk.Frame(col_frame, bg=THEME["card"])
        left_col.pack(side="left", fill="both", expand=True)
        right_col.pack(side="right", fill="both", expand=True)

        for i, (key, label, default, tip, dimensional) in enumerate(fields):
            target = left_col if i % 2 == 0 else right_col
            initial = old_values.get(key, default)
            var, _ = self._labeled_entry(target, label, initial, tip, width=12)
            var._dimensional = dimensional  # tag for unit conversion
            self.input_vars[key] = var

        self._draw_schematic()

    # --------------------------------------------------------
    def _reset_inputs(self):
        """Reset all input fields to their default values."""
        self._refresh_inputs(preserve=False)

    def _save_setup(self):
        """Write the current test configuration to a JSON file."""
        try:
            params = {k: v.get() for k, v in self.input_vars.items()}
        except Exception as e:
            messagebox.showerror("Save error", str(e), parent=self)
            return
        setup = {
            "app": "BiaxialStrengthCalculator",
            "version": 1,
            "test_method": self.test_method.get(),
            "unit_force": self.unit_force.get(),
            "unit_length": self.unit_length.get(),
            "parameters": params,
        }
        path = filedialog.asksaveasfilename(
            title="Save setup",
            defaultextension=".json",
            filetypes=[("Setup files", "*.json"), ("All files", "*.*")],
            initialfile="biaxial_setup.json")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(setup, f, indent=2)
        except Exception as e:
            messagebox.showerror("Save error", str(e), parent=self)
            return
        messagebox.showinfo("Setup saved",
                            f"Saved to {Path(path).name}", parent=self)

    def _load_setup(self):
        """Read a previously saved JSON setup file and apply it."""
        path = filedialog.askopenfilename(
            title="Load setup",
            filetypes=[("Setup files", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                setup = json.load(f)
            if setup.get("app") != "BiaxialStrengthCalculator":
                if not messagebox.askyesno(
                        "Unrecognised file",
                        "This file does not appear to be a valid setup file. "
                        "Attempt to load it anyway?", parent=self):
                    return
        except Exception as e:
            messagebox.showerror("Load error", str(e), parent=self)
            return

        # Apply method first so the right field set is built
        method = setup.get("test_method", self.test_method.get())
        self.test_method.set(method)
        self.unit_force.set(setup.get("unit_force", self.unit_force.get()))
        self.unit_length.set(setup.get("unit_length", self.unit_length.get()))

        self._refresh_inputs(preserve=False)

        # Now overwrite the freshly built defaults with the loaded values
        params = setup.get("parameters", {})
        for k, val in params.items():
            if k in self.input_vars:
                self.input_vars[k].set(str(val))

        messagebox.showinfo("Setup loaded",
                            f"Loaded {Path(path).name}", parent=self)

    # --------------------------------------------------------
    def _draw_schematic(self):
        c = self.schematic
        c.delete("all")
        c.update_idletasks()
        W = c.winfo_width() or 460
        H = c.winfo_height() or 280
        if W < 20 or H < 20:  # not yet laid out
            return

        # Reserve ~70 px on the right for the "specimen" label
        label_reserve = 70
        usable_w = max(160, W - label_reserve - 20)

        # Disc occupies most of the usable width
        disc_w = min(usable_w, 320)
        disc_h = 16
        cx = 10 + disc_w // 2  # disc centred in the usable area (left-ish)
        cy = H // 2
        disc_top = cy - disc_h // 2
        disc_bot = cy + disc_h // 2
        disc_left = cx - disc_w // 2
        disc_right = cx + disc_w // 2

        method = self.test_method.get()

        # --- Common: specimen disc (side view) ---
        c.create_rectangle(disc_left, disc_top, disc_right, disc_bot,
                           fill="#E8EEF5", outline="#7A99B8", width=1.5)
        c.create_text(disc_right + 12, cy, anchor="w",
                      text="specimen", fill=THEME["text_muted"],
                      font=(FONT_FALLBACK, 9, "italic"))

        ball_r = 7  # radius of all schematic balls
        # Scale ball separation to disc size
        sup_offset = int(disc_w * 0.30)
        load_offset = int(disc_w * 0.14)

        if method.startswith("Ring-on-Ring"):
            # Support ring positions (sit just below the disc bottom)
            sup_x = (cx - sup_offset, cx + sup_offset)
            sup_cy = disc_bot + ball_r + 1  # ball centre
            for x in sup_x:
                c.create_oval(x - ball_r, sup_cy - ball_r,
                              x + ball_r, sup_cy + ball_r,
                              fill=THEME["accent"], outline="")
            # D_S label centred under the support ring with leader lines
            ds_label_y = sup_cy + ball_r + 22
            c.create_line(sup_x[0], sup_cy + ball_r + 2,
                          sup_x[0], ds_label_y - 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_line(sup_x[1], sup_cy + ball_r + 2,
                          sup_x[1], ds_label_y - 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_line(sup_x[0], ds_label_y - 6,
                          sup_x[1], ds_label_y - 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_text(cx, ds_label_y, text="D_S",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Load ring positions (sit just above the disc top)
            load_x = (cx - load_offset, cx + load_offset)
            load_cy = disc_top - ball_r - 1
            for x in load_x:
                c.create_oval(x - ball_r, load_cy - ball_r,
                              x + ball_r, load_cy + ball_r,
                              fill=THEME["accent_press"], outline="")
            # D_L label above with leader lines
            dl_label_y = load_cy - ball_r - 22
            c.create_line(load_x[0], load_cy - ball_r - 2,
                          load_x[0], dl_label_y + 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_line(load_x[1], load_cy - ball_r - 2,
                          load_x[1], dl_label_y + 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_line(load_x[0], dl_label_y + 6,
                          load_x[1], dl_label_y + 6,
                          fill=THEME["text_subtle"], width=1)
            c.create_text(cx, dl_label_y, text="D_L",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Force arrow — well above the D_L bracket so nothing overlaps
            arrow_top = dl_label_y - 38
            arrow_bot = dl_label_y - 10
            c.create_line(cx, arrow_top, cx, arrow_bot,
                          arrow="last", fill=THEME["text"], width=2)
            c.create_text(cx + 14, (arrow_top + arrow_bot) // 2,
                          text="F", fill=THEME["text"],
                          font=(FONT_FALLBACK, 12, "bold"))

        elif method.startswith("Piston"):
            # Three support balls (bottom — show 2 in side view)
            sup_x = (cx - sup_offset, cx + sup_offset)
            sup_cy = disc_bot + ball_r + 1
            for x in sup_x:
                c.create_oval(x - ball_r, sup_cy - ball_r,
                              x + ball_r, sup_cy + ball_r,
                              fill=THEME["accent"], outline="")
            c.create_text(cx, sup_cy + ball_r + 16,
                          text="3 support balls (r1)",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Piston: long cylindrical pin with flat tip sitting on the disc
            piston_w = 10      # pin diameter in the side view
            piston_h = 55      # pin length
            p_bot = disc_top   # flat tip touches disc top
            p_top = p_bot - piston_h
            # Main cylindrical body
            c.create_rectangle(cx - piston_w // 2, p_top,
                               cx + piston_w // 2, p_bot,
                               fill="#B0B0B0", outline="#555", width=1)
            # Subtle vertical highlight for cylindrical shading
            c.create_line(cx - piston_w // 2 + 2, p_top + 2,
                          cx - piston_w // 2 + 2, p_bot - 2,
                          fill="#D8D8D8", width=1)
            # Label to the right, joined by a thin leader to the tip
            piston_label_x = cx + piston_w // 2 + 38
            piston_label_y = p_bot - 8
            c.create_line(cx + piston_w // 2 + 2, piston_label_y,
                          piston_label_x - 4, piston_label_y,
                          fill=THEME["text_subtle"], width=1)
            c.create_text(piston_label_x, piston_label_y,
                          text="piston (r2)", anchor="w",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Force arrow above the piston
            arrow_top = p_top - 38
            arrow_bot = p_top - 4
            c.create_line(cx, arrow_top, cx, arrow_bot,
                          arrow="last", fill=THEME["text"], width=2)
            c.create_text(cx + 14, (arrow_top + arrow_bot) // 2,
                          text="F", fill=THEME["text"],
                          font=(FONT_FALLBACK, 12, "bold"))

        else:  # B3B
            # Three support balls (bottom — show 2 in side view)
            sup_x = (cx - sup_offset, cx + sup_offset)
            sup_cy = disc_bot + ball_r + 1
            for x in sup_x:
                c.create_oval(x - ball_r, sup_cy - ball_r,
                              x + ball_r, sup_cy + ball_r,
                              fill=THEME["accent"], outline="")
            c.create_text(cx, sup_cy + ball_r + 16,
                          text="3 support balls",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Loading ball on top, centred, touching the disc top
            load_cy = disc_top - ball_r - 1
            c.create_oval(cx - ball_r, load_cy - ball_r,
                          cx + ball_r, load_cy + ball_r,
                          fill=THEME["accent_press"], outline="")
            # Label to the right with a leader line
            load_label_x = cx + ball_r + 38
            c.create_line(cx + ball_r + 2, load_cy,
                          load_label_x - 4, load_cy,
                          fill=THEME["text_subtle"], width=1)
            c.create_text(load_label_x, load_cy,
                          text="loading ball", anchor="w",
                          fill=THEME["text_muted"],
                          font=(FONT_FALLBACK, 9))

            # Force arrow well above the loading ball
            arrow_top = load_cy - ball_r - 40
            arrow_bot = load_cy - ball_r - 6
            c.create_line(cx, arrow_top, cx, arrow_bot,
                          arrow="last", fill=THEME["text"], width=2)
            c.create_text(cx + 14, (arrow_top + arrow_bot) // 2,
                          text="F", fill=THEME["text"],
                          font=(FONT_FALLBACK, 12, "bold"))

    # --------------------------------------------------------
    def _get_inputs_SI(self):
        """Read entries, convert to SI, return dict and the raw force/lengths."""
        force_unit = self.unit_force.get()
        length_unit = self.unit_length.get()
        f_scale = 1000.0 if force_unit == "kN" else 1.0
        l_scale = 1e-3 if length_unit == "mm" else 1.0

        out = {}
        for key, var in self.input_vars.items():
            try:
                v = float(var.get())
            except ValueError:
                raise ValueError(f"Invalid number for '{key}'.")
            if key == "F":
                out[key] = v * f_scale
            elif key in ("nu", "f"):
                out[key] = v
            else:
                out[key] = v * l_scale
        return out

    def _calculate(self):
        try:
            vals = self._get_inputs_SI()
            method = self.test_method.get()
            if method.startswith("Ring-on-Ring"):
                sigma = StrengthCalculator.ring_on_ring(
                    vals["F"], vals["h"], vals["nu"],
                    vals["DS"], vals["DL"], vals["D"])
            elif method.startswith("Piston"):
                sigma = StrengthCalculator.piston_on_three_balls(
                    vals["F"], vals["h"], vals["nu"],
                    vals["r1"], vals["r2"], vals["D"] / 2.0)
            else:
                sigma = StrengthCalculator.ball_on_three_balls(
                    vals["F"], vals["h"], vals["f"])
            sigma_MPa = sigma / 1e6
            self.result_value.config(text=f"{sigma_MPa:,.2f}")
            self.result_detail.config(
                text=(f"Method: {method}\n"
                      f"Force at failure: {vals['F']:.3f} N    "
                      f"Thickness: {vals['h']*1000:.3f} mm"))
            self._last_sigma_MPa = sigma_MPa
            self._last_force_N = vals["F"]
            return sigma_MPa
        except Exception as e:
            messagebox.showerror("Calculation error", str(e), parent=self)
            return None

    def _add_to_batch(self):
        sigma_MPa = self._calculate()
        if sigma_MPa is None:
            return
        self.specimens.append({
            "force_N": self._last_force_N,
            "sigma_MPa": sigma_MPa,
            "method": self.test_method.get(),
            "valid": True,
        })
        self._refresh_batch_table()
        messagebox.showinfo("Added",
                            f"Specimen added to batch.\nσ = {sigma_MPa:.2f} MPa\n"
                            f"Total specimens: {len(self.specimens)}",
                            parent=self)

    # --------------------------------------------------------
    def _build_weibull_tab(self):
        container = tk.Frame(self.tab_weibull, bg=THEME["bg"])
        container.pack(fill="both", expand=True)

        # Left: data + buttons
        left = tk.Frame(container, bg=THEME["bg"])
        left.pack(side="left", fill="both", expand=False, padx=(0, 10))

        data_outer, data_inner = self._card(left, "Specimen batch")
        data_outer.pack(fill="both", expand=True, pady=(0, 12))

        # Toolbar
        tb = tk.Frame(data_inner, bg=THEME["card"])
        tb.pack(fill="x", padx=18, pady=(0, 8))
        WinButton(tb, "Import CSV", command=self._import_csv,
                  style="subtle", width=120, bg=THEME["card"]).pack(side="left")
        WinButton(tb, "Add manually", command=self._add_manual_dialog,
                  style="subtle", width=120, bg=THEME["card"]).pack(side="left", padx=6)
        WinButton(tb, "Remove", command=self._remove_selected,
                  style="subtle", width=90, bg=THEME["card"]).pack(side="left")
        WinButton(tb, "Clear all", command=self._clear_batch,
                  style="subtle", width=90, bg=THEME["card"]).pack(side="left", padx=6)

        # Treeview
        tv_frame = tk.Frame(data_inner, bg=THEME["card"])
        tv_frame.pack(fill="both", expand=True, padx=18, pady=(0, 14))
        cols = ("idx", "force", "sigma", "method")
        self.tree = ttk.Treeview(tv_frame, columns=cols, show="headings",
                                 style="Win.Treeview", height=14)
        self.tree.heading("idx", text="#")
        self.tree.heading("force", text="F (N)")
        self.tree.heading("sigma", text="σ (MPa)")
        self.tree.heading("method", text="Method")
        self.tree.column("idx", width=40, anchor="center")
        self.tree.column("force", width=90, anchor="e")
        self.tree.column("sigma", width=90, anchor="e")
        self.tree.column("method", width=240, anchor="w")
        self.tree.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(tv_frame, orient="vertical", command=self.tree.yview)
        sb.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=sb.set)

        # Right: results & plot
        right = tk.Frame(container, bg=THEME["bg"])
        right.pack(side="right", fill="both", expand=True, padx=(10, 0))

        stats_outer, stats_inner = self._card(right, "Weibull parameters")
        stats_outer.pack(fill="x", pady=(0, 12))

        self.stats_labels = {}
        for key, label, tip in [
            ("n", "Sample size, n",
             "Number of valid specimens used in the fit."),
            ("m", "Weibull modulus, m",
             "Shape parameter — higher m means more uniform strength "
             "(narrower distribution). Typical ceramics: 5–20."),
            ("s0", "Characteristic strength, σ₀ (MPa)",
             "Scale parameter — the strength at which 63.2% of specimens "
             "would have failed."),
            ("smean", "Mean strength (MPa)", "Arithmetic mean of the dataset."),
            ("sstd", "Std. deviation (MPa)", "Sample standard deviation."),
        ]:
            row = tk.Frame(stats_inner, bg=THEME["card"])
            row.pack(fill="x", padx=18, pady=3)
            lbl = tk.Label(row, text=label, bg=THEME["card"],
                           fg=THEME["text_muted"],
                           font=(FONT_FALLBACK, 10), anchor="w")
            lbl.pack(side="left")
            Tooltip(lbl, tip)
            val = tk.Label(row, text="—", bg=THEME["card"],
                           fg=THEME["text"],
                           font=(FONT_FALLBACK, 11, "bold"))
            val.pack(side="right")
            self.stats_labels[key] = val

        btn_row = tk.Frame(right, bg=THEME["bg"])
        btn_row.pack(fill="x", pady=(0, 12))
        WinButton(btn_row, "Run Weibull fit", command=self._run_weibull,
                  style="accent", width=160, bg=THEME["bg"]).pack(side="left")
        WinButton(btn_row, "Export results", command=self._export_results,
                  style="subtle", width=140, bg=THEME["bg"]).pack(side="left", padx=8)
        WinButton(btn_row, "Export plot data", command=self._export_plot_data,
                  style="subtle", width=140, bg=THEME["bg"]).pack(side="left")

        plot_outer, plot_inner = self._card(right, "Weibull plot")
        plot_outer.pack(fill="both", expand=True)
        self.plot_frame = tk.Frame(plot_inner, bg=THEME["card"])
        self.plot_frame.pack(fill="both", expand=True, padx=18, pady=(4, 14))

        if not HAS_PLOTTING:
            tk.Label(self.plot_frame,
                     text="numpy and matplotlib are required for Weibull fitting and plotting.",
                     bg=THEME["card"], fg=THEME["error"],
                     font=(FONT_FALLBACK, 10, "italic")).pack(pady=20)

    # --------------------------------------------------------
    def _refresh_batch_table(self):
        for it in self.tree.get_children():
            self.tree.delete(it)
        for i, s in enumerate(self.specimens, 1):
            self.tree.insert("", "end", values=(
                i,
                f"{s['force_N']:.2f}",
                f"{s['sigma_MPa']:.2f}",
                s["method"],
            ))

    def _remove_selected(self):
        sel = self.tree.selection()
        if not sel:
            return
        idxs = sorted([self.tree.index(item) for item in sel], reverse=True)
        for i in idxs:
            del self.specimens[i]
        self._refresh_batch_table()

    def _clear_batch(self):
        if not self.specimens:
            return
        if messagebox.askyesno("Clear all",
                               "Remove all specimens from the batch?",
                               parent=self):
            self.specimens.clear()
            self._refresh_batch_table()

    def _add_manual_dialog(self):
        dlg = tk.Toplevel(self)
        dlg.title("Add specimen")
        dlg.configure(bg=THEME["card"])
        dlg.transient(self)
        dlg.grab_set()
        dlg.geometry("340x180")
        tk.Label(dlg, text="Add specimen manually",
                 bg=THEME["card"], fg=THEME["text"],
                 font=(FONT_FALLBACK, 12, "bold")).pack(anchor="w",
                                                        padx=18, pady=(14, 8))
        v_sigma = tk.StringVar()
        v_force = tk.StringVar(value="0")
        for label, var in [("Strength σ (MPa)", v_sigma),
                           ("Failure load F (N) [optional]", v_force)]:
            row = tk.Frame(dlg, bg=THEME["card"])
            row.pack(fill="x", padx=18, pady=4)
            tk.Label(row, text=label, bg=THEME["card"], fg=THEME["text"],
                     font=(FONT_FALLBACK, 10)).pack(side="left")
            tk.Entry(row, textvariable=var, bg=THEME["input_bg"],
                     relief="solid", bd=1,
                     highlightthickness=1,
                     highlightbackground=THEME["input_border"],
                     highlightcolor=THEME["input_focus"]
                     ).pack(side="right", ipady=3)

        def add():
            try:
                s = float(v_sigma.get())
                f = float(v_force.get() or 0)
            except ValueError:
                messagebox.showerror("Invalid", "Enter numeric values.", parent=dlg)
                return
            self.specimens.append({"force_N": f, "sigma_MPa": s,
                                   "method": "Manual entry", "valid": True})
            self._refresh_batch_table()
            dlg.destroy()

        btns = tk.Frame(dlg, bg=THEME["card"])
        btns.pack(fill="x", padx=18, pady=12)
        WinButton(btns, "Add", command=add,
                  style="accent", width=100, bg=THEME["card"]).pack(side="right")
        WinButton(btns, "Cancel", command=dlg.destroy,
                  style="subtle", width=100, bg=THEME["card"]
                  ).pack(side="right", padx=8)

    def _import_csv(self):
        path = filedialog.askopenfilename(
            title="Import strength data",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        added = 0
        try:
            with open(path, newline="", encoding="utf-8-sig") as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row:
                        continue
                    try:
                        sigma = float(row[0])
                    except ValueError:
                        continue  # skip headers
                    force = 0.0
                    if len(row) > 1:
                        try:
                            force = float(row[1])
                        except ValueError:
                            force = 0.0
                    self.specimens.append({"force_N": force,
                                           "sigma_MPa": sigma,
                                           "method": "Imported",
                                           "valid": True})
                    added += 1
        except Exception as e:
            messagebox.showerror("Import error", str(e), parent=self)
            return
        self._refresh_batch_table()
        messagebox.showinfo("Import", f"Imported {added} specimens.", parent=self)

    # --------------------------------------------------------
    def _run_weibull(self):
        if not HAS_PLOTTING:
            messagebox.showerror("Missing libraries",
                                 "numpy and matplotlib are required for this feature.",
                                 parent=self)
            return
        strengths = [s["sigma_MPa"] for s in self.specimens if s["valid"]]
        if len(strengths) < 2:
            messagebox.showwarning("Not enough data",
                                   "Add at least 2 specimens.", parent=self)
            return
        try:
            m, s0, n = WeibullAnalysis.fit(strengths)
        except Exception as e:
            messagebox.showerror("Fit failed", str(e), parent=self)
            return

        arr = np.asarray(strengths)
        self.stats_labels["n"].config(text=f"{n}")
        self.stats_labels["m"].config(text=f"{m:.2f}")
        self.stats_labels["s0"].config(text=f"{s0:.2f}")
        self.stats_labels["smean"].config(text=f"{arr.mean():.2f}")
        self.stats_labels["sstd"].config(text=f"{arr.std(ddof=1):.2f}")

        # Store fit results for plot-data export
        self._weibull_m = m
        self._weibull_s0 = s0

        self._draw_weibull_plot(arr, m, s0)

    def _draw_weibull_plot(self, strengths, m, s0):
        for w in self.plot_frame.winfo_children():
            w.destroy()
        fig = Figure(figsize=(5.4, 3.8), dpi=100, facecolor=THEME["card"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(THEME["card"])

        x_sorted, P = WeibullAnalysis.probability_of_failure(strengths)
        # Linearised Weibull: ln(ln(1/(1-P))) vs ln(sigma)
        ln_sigma = np.log(x_sorted)
        ln_lnP = np.log(np.log(1.0 / (1.0 - P)))

        # Store the plot data for export
        self._weibull_plot_data = {
            "sigma_sorted": x_sorted,
            "P": P,
            "ln_sigma": ln_sigma,
            "ln_lnP": ln_lnP,
        }

        ax.scatter(ln_sigma, ln_lnP, color=THEME["accent"],
                   s=40, edgecolor="white", linewidth=1, zorder=3,
                   label="Specimens")

        # Fit line from MLE m, s0
        x_line = np.linspace(ln_sigma.min() - 0.05, ln_sigma.max() + 0.05, 100)
        y_line = m * x_line - m * np.log(s0)
        ax.plot(x_line, y_line, color=THEME["accent_press"],
                linewidth=2, label=f"MLE fit (m = {m:.2f})")

        ax.set_xlabel("ln(σ)  [σ in MPa]", fontsize=9)
        ax.set_ylabel("ln(ln(1/(1−P)))", fontsize=9)
        ax.set_title("Weibull plot", fontsize=10, color=THEME["text"])
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
        ax.tick_params(labelsize=8)
        for spine in ax.spines.values():
            spine.set_color(THEME["card_border"])
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    def _export_results(self):
        if not self.specimens:
            messagebox.showinfo("Nothing to export", "Batch is empty.", parent=self)
            return
        path = filedialog.asksaveasfilename(
            title="Export results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["#", "Force (N)", "Strength (MPa)", "Method"])
            for i, s in enumerate(self.specimens, 1):
                w.writerow([i, f"{s['force_N']:.3f}",
                            f"{s['sigma_MPa']:.3f}", s["method"]])
            # Stats
            for key in ("n", "m", "s0", "smean", "sstd"):
                w.writerow([])
                break
            w.writerow([])
            w.writerow(["Weibull statistics"])
            for key, label in [("n", "n"), ("m", "Weibull modulus m"),
                               ("s0", "Characteristic strength sigma_0 (MPa)"),
                               ("smean", "Mean strength (MPa)"),
                               ("sstd", "Std deviation (MPa)")]:
                w.writerow([label, self.stats_labels[key].cget("text")])
        messagebox.showinfo("Export", f"Saved to {Path(path).name}", parent=self)

    def _export_plot_data(self):
        """Export Weibull plot coordinates and fit line to CSV."""
        if not hasattr(self, "_weibull_plot_data") or self._weibull_plot_data is None:
            messagebox.showinfo("No plot data",
                                "Run a Weibull fit first.", parent=self)
            return
        path = filedialog.asksaveasfilename(
            title="Export Weibull plot data",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="weibull_plot_data.csv")
        if not path:
            return

        d = self._weibull_plot_data
        m = self._weibull_m
        s0 = self._weibull_s0
        n = len(d["sigma_sorted"])

        # Build the fit line over the same ln(sigma) range
        ln_sig = d["ln_sigma"]
        x_fit = np.linspace(float(ln_sig.min()) - 0.05,
                            float(ln_sig.max()) + 0.05, 100)
        y_fit = m * x_fit - m * np.log(s0)

        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            # --- Header block: metadata for whoever opens this ---
            w.writerow(["# Weibull plot data exported from Biaxial Flexural "
                        "Strength Calculator v1.0"])
            w.writerow([f"# Weibull modulus m = {m:.4f}"])
            w.writerow([f"# Characteristic strength sigma_0 = {s0:.4f} MPa"])
            w.writerow([f"# Sample size n = {n}"])
            w.writerow([f"# Fit line equation: y = m * x - m * ln(sigma_0)  "
                        f"i.e.  y = {m:.4f} * x - {m * np.log(s0):.4f}"])
            w.writerow([])

            # --- Section 1: Specimen data points ---
            w.writerow(["# SECTION 1: Specimen data points (for scatter plot)"])
            w.writerow(["# X-axis: ln(sigma)  [sigma in MPa]"])
            w.writerow(["# Y-axis: ln(ln(1/(1-P)))  [dimensionless]"])
            w.writerow([])
            w.writerow(["Rank_i",
                        "Strength_sigma (MPa)",
                        "Probability_of_failure_P",
                        "ln(sigma)",
                        "ln(ln(1/(1-P)))"])
            for i in range(n):
                w.writerow([i + 1,
                            f"{float(d['sigma_sorted'][i]):.4f}",
                            f"{float(d['P'][i]):.6f}",
                            f"{float(d['ln_sigma'][i]):.6f}",
                            f"{float(d['ln_lnP'][i]):.6f}"])
            w.writerow([])

            # --- Section 2: MLE fit line ---
            w.writerow(["# SECTION 2: MLE fit line (for line plot overlay)"])
            w.writerow(["# X-axis: ln(sigma)  [sigma in MPa]"])
            w.writerow(["# Y-axis: ln(ln(1/(1-P)))  [dimensionless]"])
            w.writerow([])
            w.writerow(["ln(sigma)_fit", "ln(ln(1/(1-P)))_fit"])
            for i in range(len(x_fit)):
                w.writerow([f"{float(x_fit[i]):.6f}",
                            f"{float(y_fit[i]):.6f}"])

        messagebox.showinfo("Export",
                            f"Plot data saved to {Path(path).name}\n\n"
                            "Contents:\n"
                            "  1. Specimen data points\n"
                            "  2. MLE fit line (100 points)\n\n"
                            "Axis titles and units are included as "
                            "comment headers.",
                            parent=self)

    # --------------------------------------------------------
    def _build_about_tab(self):
        outer, inner = self._card(self.tab_about, "About")
        outer.pack(fill="both", expand=True, padx=0, pady=0)
        text = (
            "Biaxial Flexural Strength Calculator  ·  v1.0\n\n"
            "Author: Dr Thanos Goulas\n"
            "Contact: thanosgoulas@outlook.com\n\n"
            "Feedback, bug reports and suggestions are welcome. "
            "Please open an issue on the GitHub repository or get "
            "in touch by email.\n\n"
            "This tool calculates the biaxial flexural strength of "
            "ceramic disc specimens using the closed-form stress "
            "solutions for three standard test configurations:\n\n"
            "  •  Ring-on-Ring (ASTM C1499)\n"
            "  •  Piston-on-three-balls (ISO 6872)\n"
            "  •  Ball-on-three-balls (B3B)\n\n"
            "A two-parameter Weibull statistical module is included "
            "for batch analysis of strength data. The Weibull modulus "
            "and characteristic strength are fitted by maximum "
            "likelihood estimation. Specimen data, Weibull parameters "
            "and plot coordinates can be exported to CSV for further "
            "processing in external plotting software.\n\n"
            "Stress output is in MPa. All calculations are performed "
            "internally in SI units. The exact equations used are "
            "documented in the 'Formulas & references' tab.\n\n"
            "─────────────────────────────────────────────────────────\n"
            "Copyright (c) 2026 Dr Thanos Goulas.\n"
            "Released under the MIT License.\n\n"
            "This software is provided as an engineering aid. The "
            "accuracy of calculated values depends on the validity "
            "of the input parameters and the underlying assumptions "
            "of the chosen test standard. The user is responsible "
            "for verifying results against the relevant standard and "
            "for confirming the appropriateness of the test method "
            "for their specific application.\n\n"
            "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY "
            "OF ANY KIND, EXPRESS OR IMPLIED."
        )
        tk.Label(inner, text=text, bg=THEME["card"], fg=THEME["text"],
                 font=(FONT_FALLBACK, 10), justify="left",
                 wraplength=900).pack(anchor="w", padx=20, pady=20)

    # --------------------------------------------------------
    def _build_formulas_tab(self):
        """Populate the Formulas & references tab with the implemented equations."""
        outer, inner = self._card(self.tab_formulas, "Implemented formulas")
        outer.pack(fill="both", expand=True, padx=0, pady=0)

        # Scrollable container so longer content stays accessible
        canvas = tk.Canvas(inner, bg=THEME["card"], highlightthickness=0)
        sb = ttk.Scrollbar(inner, orient="vertical", command=canvas.yview)
        body = tk.Frame(canvas, bg=THEME["card"])
        body.bind("<Configure>",
                  lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=body, anchor="nw")
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side="left", fill="both", expand=True, padx=(18, 0), pady=(4, 14))
        sb.pack(side="right", fill="y", pady=(4, 14))

        def section(title):
            tk.Label(body, text=title, bg=THEME["card"], fg=THEME["accent"],
                     font=(FONT_FALLBACK, 12, "bold")).pack(anchor="w", pady=(14, 4))

        def para(text, mono=False):
            font = ("Consolas", 10) if mono else (FONT_FALLBACK, 10)
            tk.Label(body, text=text, bg=THEME["card"], fg=THEME["text"],
                     font=font, justify="left", wraplength=860
                     ).pack(anchor="w", pady=2)

        def muted(text):
            tk.Label(body, text=text, bg=THEME["card"], fg=THEME["text_muted"],
                     font=(FONT_FALLBACK, 9, "italic"), justify="left",
                     wraplength=860).pack(anchor="w", pady=(2, 6))

        # --- Common ---
        section("Symbols and units")
        para("F   = force at failure  [N]\n"
             "h   = specimen thickness  [m]\n"
             "ν   = Poisson's ratio  [–]\n"
             "σ   = maximum tensile stress on the lower (tensile) face  [Pa]\n"
             "Lengths are entered in mm or m (selectable). Force is entered "
             "in N or kN (selectable). All formulas below are in SI units; "
             "inputs are converted internally and σ is reported in MPa.",
             mono=False)

        # --- Ring on Ring ---
        section("1. Ring-on-Ring  (ASTM C1499)")
        para("σ = (3·F / (2·π·h²)) · [ (1−ν)·(D_S² − D_L²) / (2·D²)  +  (1+ν)·ln(D_S / D_L) ]",
             mono=True)
        para("D_S = support ring diameter,  D_L = load ring diameter,  "
             "D = specimen diameter.")
        muted("Reference: ASTM C1499 — Standard Test Method for Monotonic "
              "Equibiaxial Flexural Strength of Advanced Ceramics at Ambient "
              "Temperature. Stress is the maximum equibiaxial tensile stress "
              "on the lower face within the load ring.")

        # --- Piston on three balls ---
        section("2. Piston-on-three-balls  (ISO 6872)")
        para("σ = −0.2387 · F · (X − Y) / h²", mono=True)
        para("with", mono=False)
        para("X = (1 + ν) · ln((r₂ / r₃)²)  +  (1 − ν) · (r₂² − r₃²) / (2·r₁²)\n"
             "Y = (1 + ν) · [1 + ln((r₁ / r₃)²)]  +  (1 − ν) · (r₁² − r₃²) / r₃²",
             mono=True)
        para("r₁ = radius of the support circle through the centres of the "
             "three balls\n"
             "r₂ = contact radius of the piston tip on the disc\n"
             "r₃ = specimen radius (= D / 2)")
        muted("Reference: ISO 6872 — Dentistry — Ceramic materials. The form "
              "above is the ISO 6872 / Wachtman expression for monolayered "
              "discs. ISO 6872 recommends r₂ = h/3 when the piston contact "
              "radius is not known.")

        # --- B3B ---
        section("3. Ball-on-three-balls  (B3B)")
        para("σ = f · F / h²", mono=True)
        para("f is a dimensionless geometry factor obtained from FEA. It "
             "depends on the specimen radius, the support ball radius and "
             "Poisson's ratio. Typical values lie in the range 1.2 – 2.0.")
        muted("References: Börger A., Supancic P., Danzer R. The ball on "
              "three balls test for strength testing of brittle discs: "
              "stress distribution in the disc. J. Eur. Ceram. Soc. "
              "22 (2002) 1425–1436.  And: Börger A., Supancic P., Danzer R. "
              "The ball on three balls test for strength testing of brittle "
              "discs: Part II — analysis of possible errors in the strength "
              "determination. J. Eur. Ceram. Soc. 24 (2004) 2917–2928.")

        # --- Weibull ---
        section("4. Two-parameter Weibull statistics")
        para("Probability of failure (cumulative distribution):\n"
             "P(σ) = 1 − exp[ −(σ / σ₀)^m ]", mono=True)
        para("m = Weibull modulus (shape parameter)\n"
             "σ₀ = characteristic strength (the stress at which 63.2 % of "
             "specimens have failed)")
        para("m and σ₀ are determined by maximum likelihood estimation (MLE). "
             "m is obtained by Newton–Raphson iteration on the MLE score "
             "equation:", mono=False)
        para("Σ σᵢ^m · ln(σᵢ) / Σ σᵢ^m  −  (1/n)·Σ ln(σᵢ)  −  1/m  =  0",
             mono=True)
        para("σ₀ then follows from:", mono=False)
        para("σ₀ = ( (1/n) · Σ σᵢ^m )^(1/m)", mono=True)
        para("On the Weibull plot, the empirical probability of failure "
             "uses the median rank estimator:", mono=False)
        para("Pᵢ = (i − 0.5) / n", mono=True)
        para("for the i-th specimen (i = 1 … n) when strengths are sorted "
             "in ascending order. The linearised plot displays ln(ln(1/(1−P))) "
             "versus ln(σ), with the MLE fit line ln(ln(1/(1−P))) = m·ln(σ) "
             "− m·ln(σ₀) overlaid.")
        muted("References: Weibull W. A statistical distribution function of "
              "wide applicability. J. Appl. Mech. 18 (1951) 293–297.  See also "
              "ENV 843-5 / EN 843-5 for guidance on Weibull analysis of "
              "ceramic strength data, including recommended sample sizes.")


# ============================================================
if __name__ == "__main__":
    app = BiaxialApp()
    app.mainloop()