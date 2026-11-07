"""
SeqAlign+ : Kivy Desktop Sequence Alignment Application
======================================================

A cross-platform graphical desktop application built with Kivy for interactive
pairwise biological sequence alignment.

Features:
- Global Alignment (Needleman-Wunsch) with interactive dynamic programming matrix grid.
- Local Alignment (Smith-Waterman) with custom scoring scheme inputs.
- Multi-screen navigation with fluid slide transitions.
- Visual display of scoring parameters (Match, Mismatch, Gap penalty).

Author: Azheem Syahmee
License: MIT
"""

import sys
import os

# Ensure package root is in path when executed directly
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, "../.."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

try:
    from kivy.app import App
    from kivy.uix.screenmanager import ScreenManager, Screen, SlideTransition
    from kivy.uix.button import Button
    from kivy.uix.label import Label
    from kivy.uix.boxlayout import BoxLayout
    from kivy.uix.gridlayout import GridLayout
    from kivy.uix.textinput import TextInput
    from kivy.uix.scrollview import ScrollView
except ImportError:
    print("Warning: Kivy is not installed. To run this GUI, install via: pip install kivy")

from src.algorithms.needleman_wunsch import needleman_wunsch
from src.algorithms.smith_waterman import smith_waterman


class HomePage(Screen):
    """Initial landing page with navigation options."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        layout = BoxLayout(orientation="vertical", padding=20, spacing=15)

        welcome_label = Label(
            text="SeqAlign+\nBiological Sequence Alignment Suite",
            font_size=24,
            halign="center",
            size_hint_y=None,
            height=80,
        )

        desc_label = Label(
            text="Select an alignment mode below to begin pairwise sequence analysis:",
            font_size=14,
            halign="center",
            size_hint_y=None,
            height=40,
        )

        nw_btn = Button(
            text="Needleman-Wunsch (Global Alignment)",
            size_hint_y=None,
            height=50,
            background_color=(0.2, 0.6, 0.8, 1),
        )
        nw_btn.bind(on_press=self.go_to_needleman)

        sw_btn = Button(
            text="Smith-Waterman (Local Alignment)",
            size_hint_y=None,
            height=50,
            background_color=(0.2, 0.7, 0.4, 1),
        )
        sw_btn.bind(on_press=self.go_to_waterman)

        info_btn = Button(
            text="About & Algorithms Info",
            size_hint_y=None,
            height=45,
            background_color=(0.5, 0.5, 0.5, 1),
        )
        info_btn.bind(on_press=self.go_to_settings)

        layout.add_widget(welcome_label)
        layout.add_widget(desc_label)
        layout.add_widget(nw_btn)
        layout.add_widget(sw_btn)
        layout.add_widget(info_btn)
        self.add_widget(layout)

    def go_to_needleman(self, instance):
        self.manager.current = "needleman"

    def go_to_waterman(self, instance):
        self.manager.current = "waterman"

    def go_to_settings(self, instance):
        self.manager.current = "settings"


class GlobalAlignmentPage(Screen):
    """Needleman-Wunsch global alignment page with DP matrix grid."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        layout = BoxLayout(orientation="vertical", padding=15, spacing=10)

        header = Label(
            text="Needleman-Wunsch Global Alignment",
            font_size=20,
            size_hint_y=None,
            height=40,
        )

        self.seq1_input = TextInput(
            hint_text="Enter Sequence 1 (e.g. FKHMEDPLE)",
            multiline=False,
            size_hint_y=None,
            height=40,
        )
        self.seq2_input = TextInput(
            hint_text="Enter Sequence 2 (e.g. FMDTPLNE)",
            multiline=False,
            size_hint_y=None,
            height=40,
        )

        # Scoring parameters row
        param_layout = BoxLayout(orientation="horizontal", spacing=10, size_hint_y=None, height=40)
        self.match_input = TextInput(text="1", hint_text="Match (1)", multiline=False)
        self.mismatch_input = TextInput(text="-1", hint_text="Mismatch (-1)", multiline=False)
        self.gap_input = TextInput(text="-2", hint_text="Gap (-2)", multiline=False)

        param_layout.add_widget(Label(text="Match:"))
        param_layout.add_widget(self.match_input)
        param_layout.add_widget(Label(text="Mismatch:"))
        param_layout.add_widget(self.mismatch_input)
        param_layout.add_widget(Label(text="Gap:"))
        param_layout.add_widget(self.gap_input)

        calc_btn = Button(
            text="Calculate Global Alignment",
            size_hint_y=None,
            height=45,
            background_color=(0.2, 0.6, 0.8, 1),
        )
        calc_btn.bind(on_press=self.calculate)

        self.result_label = Label(
            text="Alignment result will appear here",
            size_hint_y=None,
            height=80,
            halign="center",
        )

        # Scrollable layout for DP matrix
        self.scroll_view = ScrollView(size_hint=(1, 1))
        self.matrix_layout = GridLayout(cols=1, size_hint_y=None)
        self.matrix_layout.bind(minimum_height=self.matrix_layout.setter("height"))
        self.scroll_view.add_widget(self.matrix_layout)

        back_btn = Button(
            text="Back to Home",
            size_hint_y=None,
            height=40,
        )
        back_btn.bind(on_press=self.go_home)

        layout.add_widget(header)
        layout.add_widget(self.seq1_input)
        layout.add_widget(self.seq2_input)
        layout.add_widget(param_layout)
        layout.add_widget(calc_btn)
        layout.add_widget(self.result_label)
        layout.add_widget(self.scroll_view)
        layout.add_widget(back_btn)
        self.add_widget(layout)

    def calculate(self, instance):
        seq1 = self.seq1_input.text.strip().upper()
        seq2 = self.seq2_input.text.strip().upper()

        try:
            m = int(self.match_input.text)
            mm = int(self.mismatch_input.text)
            gap = int(self.gap_input.text)
        except ValueError:
            self.result_label.text = "Error: Match, mismatch, and gap scores must be integers."
            return

        if not seq1 or not seq2:
            self.result_label.text = "Error: Please enter both sequences."
            return

        res = needleman_wunsch(seq1, seq2, match_score=m, mismatch_penalty=mm, gap_penalty=gap)
        self.result_label.text = (
            f"Seq1: {res['aligned_seq1']}\n"
            f"      {res['similarity_str']}\n"
            f"Seq2: {res['aligned_seq2']}\n"
            f"Score: {res['score']} | Identity: {res['identity'] * 100:.1f}%"
        )

        # Render matrix in grid
        self.matrix_layout.clear_widgets()
        num_cols = len(seq2) + 2
        matrix_grid = GridLayout(cols=num_cols, size_hint_y=None, spacing=2, padding=5)
        matrix_grid.bind(minimum_height=matrix_grid.setter("height"))

        # Column headers
        matrix_grid.add_widget(Label(text=" ", size_hint=(None, None), width=35, height=30))
        matrix_grid.add_widget(Label(text=" ", size_hint=(None, None), width=35, height=30))
        for char in seq2:
            matrix_grid.add_widget(Label(text=char, size_hint=(None, None), width=35, height=30))

        # Row headers and values
        row_labels = [" "] + list(seq1)
        for i, row in enumerate(res["matrix"]):
            matrix_grid.add_widget(Label(text=row_labels[i], size_hint=(None, None), width=35, height=30))
            for score in row:
                matrix_grid.add_widget(Label(text=str(score), size_hint=(None, None), width=35, height=30))

        self.scroll_view.clear_widgets()
        self.scroll_view.add_widget(matrix_grid)

    def go_home(self, instance):
        self.manager.current = "home"


class LocalAlignmentPage(Screen):
    """Smith-Waterman local alignment page."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        layout = BoxLayout(orientation="vertical", padding=15, spacing=10)

        header = Label(
            text="Smith-Waterman Local Alignment",
            font_size=20,
            size_hint_y=None,
            height=40,
        )

        self.seq1_input = TextInput(
            hint_text="Enter Sequence 1 (e.g. ACACACTA)",
            multiline=False,
            size_hint_y=None,
            height=40,
        )
        self.seq2_input = TextInput(
            hint_text="Enter Sequence 2 (e.g. AGCACACA)",
            multiline=False,
            size_hint_y=None,
            height=40,
        )

        param_layout = BoxLayout(orientation="horizontal", spacing=10, size_hint_y=None, height=40)
        self.match_input = TextInput(text="3", hint_text="Match (3)", multiline=False)
        self.mismatch_input = TextInput(text="-1", hint_text="Mismatch (-1)", multiline=False)
        self.gap_input = TextInput(text="-2", hint_text="Gap (-2)", multiline=False)

        param_layout.add_widget(Label(text="Match:"))
        param_layout.add_widget(self.match_input)
        param_layout.add_widget(Label(text="Mismatch:"))
        param_layout.add_widget(self.mismatch_input)
        param_layout.add_widget(Label(text="Gap:"))
        param_layout.add_widget(self.gap_input)

        calc_btn = Button(
            text="Calculate Local Alignment",
            size_hint_y=None,
            height=45,
            background_color=(0.2, 0.7, 0.4, 1),
        )
        calc_btn.bind(on_press=self.calculate)

        self.result_label = Label(
            text="Alignment result will appear here",
            size_hint_y=None,
            height=120,
            halign="center",
        )

        back_btn = Button(
            text="Back to Home",
            size_hint_y=None,
            height=40,
        )
        back_btn.bind(on_press=self.go_home)

        layout.add_widget(header)
        layout.add_widget(self.seq1_input)
        layout.add_widget(self.seq2_input)
        layout.add_widget(param_layout)
        layout.add_widget(calc_btn)
        layout.add_widget(self.result_label)
        layout.add_widget(back_btn)
        self.add_widget(layout)

    def calculate(self, instance):
        seq1 = self.seq1_input.text.strip().upper()
        seq2 = self.seq2_input.text.strip().upper()

        try:
            m = int(self.match_input.text)
            mm = int(self.mismatch_input.text)
            gap = int(self.gap_input.text)
        except ValueError:
            self.result_label.text = "Error: Match, mismatch, and gap scores must be integers."
            return

        if not seq1 or not seq2:
            self.result_label.text = "Error: Please enter both sequences."
            return

        res = smith_waterman(seq1, seq2, match_score=m, mismatch_penalty=mm, gap_penalty=gap)
        if res["score"] == 0:
            self.result_label.text = "No positive local alignment found (Score = 0)."
        else:
            self.result_label.text = (
                f"Conserved Local Motif Span:\n"
                f"Seq1 [{res['start1']}..{res['end1']}]: {res['aligned_seq1']}\n"
                f"            {res['similarity_str']}\n"
                f"Seq2 [{res['start2']}..{res['end2']}]: {res['aligned_seq2']}\n"
                f"Max Score: {res['score']} | Identity: {res['identity'] * 100:.1f}%"
            )

    def go_home(self, instance):
        self.manager.current = "home"


class SettingsPage(Screen):
    """Information and algorithm parameter explanations."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        layout = BoxLayout(orientation="vertical", padding=20, spacing=15)

        title = Label(text="Algorithm Reference & Details", font_size=20, size_hint_y=None, height=40)
        info_text = (
            "• Needleman-Wunsch (1970):\n"
            "  Global alignment aligning sequences end-to-end.\n"
            "  Best for sequences of comparable lengths.\n\n"
            "• Smith-Waterman (1981):\n"
            "  Local alignment isolating highly conserved regions.\n"
            "  Zero-floor boundary prevents negative score propagation.\n\n"
            "• Scoring Schemes:\n"
            "  Match: Positive reward for identity.\n"
            "  Mismatch: Negative penalty for non-identical substitutions.\n"
            "  Gap: Negative penalty for insertion or deletion events."
        )
        info_label = Label(text=info_text, font_size=13, halign="left", valign="top")

        back_btn = Button(text="Back to Home", size_hint_y=None, height=45)
        back_btn.bind(on_press=self.go_home)

        layout.add_widget(title)
        layout.add_widget(info_label)
        layout.add_widget(back_btn)
        self.add_widget(layout)

    def go_home(self, instance):
        self.manager.current = "home"


class SeqAlignApp(App):
    """Main Kivy Application Manager."""

    def build(self):
        sm = ScreenManager(transition=SlideTransition(direction="left"))
        sm.add_widget(HomePage(name="home"))
        sm.add_widget(GlobalAlignmentPage(name="needleman"))
        sm.add_widget(LocalAlignmentPage(name="waterman"))
        sm.add_widget(SettingsPage(name="settings"))
        return sm


if __name__ == "__main__":
    SeqAlignApp().run()
