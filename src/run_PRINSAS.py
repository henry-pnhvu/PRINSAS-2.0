# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:41:41 2024

@author: NHUHA
"""

import sys
import time
import warnings
import numpy as np
import backend_functions as bf
import plot_formating as pf
import PyQt5.QtWidgets as QtWdgt
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore

import matplotlib.backends.backend_qt5agg as mpl_backend
import matplotlib.figure as mpl_figure


# Ignore UserWarning
warnings.filterwarnings("ignore", category=UserWarning)

class PRINSAS_App(QtWdgt.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QtGui.QIcon("neutron_bear.ico")) # set program logo
        self.move(100,100) # set defaut program position on screen
        self.restore_window_location() # restore program position prior to closing
        
        # obtain screen dpi required for scalling the UI elements
        screen = self.screen()
        if screen:
            dpi_default = 120.48
            dpi = screen.physicalDotsPerInch()
            resolution_default = [2560, 1440]
            resolution = [screen.size().width(), screen.size().height()]
            self.scale = min([dpi/dpi_default, 
                              resolution[0]/resolution_default[0],
                              resolution[1]/resolution_default[1]])
        else:
            self.scale = 1
        print('Screen dpi = {:.0f}'.format(dpi))
        print('Screen resolution = {:.0f} \u00D7 {:.0f}'.format(screen.size().width(),
                                                           screen.size().height()))
        print('UI Scale factor = {:.2f}\n'.format(self.scale))
        
        # assigning default fitting values and draw program ui
        self.init_variables()
        self.init_ui()
        
    # Function used for restore program position
    def closeEvent(self, event):
        """Override the close event to save geometry."""
        self.save_window_location()
        super().closeEvent(event)

    # Function used for restore program position
    def save_window_location(self):
        """Save the window geometry using QSettings."""
        settings = QtCore.QSettings("henry@pnhvu.com", "PRINSAS 2.0")
        settings.setValue("window location", self.saveGeometry())

    # Function used for restore program position
    def restore_window_location(self):
        """Restore the window geometry using QSettings."""
        settings = QtCore.QSettings("henry@pnhvu.com", "PRINSAS 2.0")
        window_location = settings.value("window location")
        if window_location:
            self.restoreGeometry(window_location)
            
    # Function used for scalling UI element
    def set_style_sheet_group_box(self, group_box, font_size, groupbox_margin_top,
                                  groupbox_padding_top, title_padding_left):
        group_box.setStyleSheet(
            f'''
            QGroupBox{{font: bold italic {round(self.scale*font_size)}px Arial;  
                       border: {round(self.scale*1)}px solid;
                       margin-top: {round(self.scale*groupbox_margin_top)}px; 
                       padding-top: {round(self.scale*groupbox_padding_top)}px;}}''' + 
            f'''
            QGroupBox:title{{subcontrol-origin: margin;
                             left: {round(self.scale*title_padding_left)}px;
                             padding: 0px {round(self.scale*5)}px 0px {round(self.scale*5)}px; }}''')
    
    # Initiate variable essential for the fitting procedure.
    def init_variables(self):
        # UI parameters
        self.input_label_col = 0
        self.input_field_col = 2
        self.input_tool_tip_col = 4
        self.input_confirm_col = 5
        self.section_spacing = round(15*self.scale)
        self.input_result_width = round(680*self.scale)
        self.confirm_button_width = round(150*self.scale)
        # Chosen file directory
        self.chosen_data_file_dir = ''
        self.chosen_data_folder_dir = ''
        self.chosen_save_folder_dir = ''
        # Large-Q background
        self.bkgrd = 0
        # Q min and Q max for analysis
        self.Qmin = 0
        self.Qmax = np.inf
        # Points per dec for result
        self.pts_per_dec = 10
        # Smoothing factor lambda
        self.lambda_ = 1
        # Contrast btw 2 phases (cm-2)
        self.contrast = 3e10
        # Sample bulk density (g/cc)
        self.density = 1
        # Dominant phase
        self.major_phase = 'solid'
        # Pore radius for SSA extrapolation
        self.r_SSA_extrapolate = 0.2
        # Number of points for SSA extrapolation
        self.num_pts_SSA_extrapolate = 7
        
        # SANS data loaded from data file
        self.QQ_origin = []
        self.IQ_origin = []
        self.dIQ_origin = []
        self.dIQ_data = []
        self.IQ_percent_dIQ = []
        self.dIQ_user = []
        
        # background-subtracted SANS data trimmed to a certain Q-max value
        self.QQ_trim = []
        self.IQ_trim = []
        self.dIQ_trim = []

        # PDSP fit result 
        self.result_phi = 0
        self.result_SSA = 0
        self.result_rr = []
        self.result_fr = []
        self.result_dVdr = []

    # Function drawing the main ui element 
    def init_ui(self):
        self.setWindowTitle('PRINSAS 2.0')
        # Set universal font to Arial with the size of 19px
        self.setStyleSheet(f'''QWidget{{font: Arial; font-size: {round(self.scale*19)}px;}}''')
        # Set input box height
        self.setStyleSheet(f'''QLineEdit{{height: {round(22*self.scale)}px;
                                          font-size: {round(self.scale*19)}px;}}''')

        # Main layout of the program
        central_widget = QtWdgt.QWidget()
        self.setCentralWidget(central_widget)
        self.main_layout = QtWdgt.QGridLayout()
        self.main_layout.setSpacing(round(10*self.scale))
        central_widget.setLayout(self.main_layout)

        # User input area
        user_input_result_panel = QtWdgt.QWidget()
        user_input_result_panel.setFixedWidth(self.input_result_width)    
        user_input_result_layout = QtWdgt.QVBoxLayout()
        user_input_result_layout.setContentsMargins(0, 0, 0, 0)
        user_input_result_panel.setLayout(user_input_result_layout)
        self.main_layout.addWidget(user_input_result_panel, 0, 0, 2, 1)
        
        # Plotting of result areas
        self.create_figure_windows()
        
        # Add file selection row to user input area
        file_selection_row = QtWdgt.QHBoxLayout()
        user_input_result_layout.addLayout(file_selection_row)
        self.create_file_selection(file_selection_row)
        
        # Add data input block to user input area
        input_group_box = QtWdgt.QGroupBox("Fitting Parameters")
        self.set_style_sheet_group_box(input_group_box, font_size=22, 
                                       groupbox_margin_top=15, 
                                       groupbox_padding_top=20, 
                                       title_padding_left=30)
        data_input_grid = QtWdgt.QGridLayout()
        data_input_grid.setColumnMinimumWidth(1,round(17*self.scale))
        data_input_grid.setColumnMinimumWidth(2,round(120*self.scale))
        data_input_grid.setColumnMinimumWidth(3,round(40*self.scale))
        data_input_grid.setVerticalSpacing(round(12*self.scale))
        data_input_grid.setRowMinimumHeight(5, round(24*self.scale))  # Spacing between lambda input and PDSP input
        input_group_box.setLayout(data_input_grid)
        user_input_result_layout.addWidget(input_group_box)
        user_input_result_layout.insertSpacing(1, self.section_spacing)
        
        # Add background input row to user input area
        self.create_bkgrd_Q_range_input(data_input_grid, row = 0)
        
        # Add number of points per decade input row and lambda input row
        # to user input area
        self.create_pts_per_dec_lambda_input(data_input_grid, row = 3)
        
        # Add PDSP input row to user input area
        self.create_PDSP_fit_input(data_input_grid, row = 6)
        
        # Add run fit button to user input area
        self.run_fit_button = QtWdgt.QPushButton("Fit PDSP Model!")
        self.format_PDSP_button()
        user_input_result_layout.addWidget(self.run_fit_button)
        user_input_result_layout.insertSpacing(3, self.section_spacing)

        # Add result display area
        result_group_box = QtWdgt.QGroupBox("PDSP Fit Results")
        self.set_style_sheet_group_box(result_group_box, font_size=22, 
                                       groupbox_margin_top=15, 
                                       groupbox_padding_top=20, 
                                       title_padding_left=30)
        user_input_result_layout.addWidget(result_group_box)
        user_input_result_layout.insertSpacing(5, self.section_spacing)
        self.result_grid = QtWdgt.QGridLayout()
        self.result_grid.setColumnMinimumWidth(1,round(5*self.scale))
        self.result_grid.setVerticalSpacing(round(12*self.scale))
        result_group_box.setLayout(self.result_grid)
        self.create_result(self.result_grid)
        user_input_result_layout.insertSpacing(7, 
                                               round(self.section_spacing*0.8))

        # Add plot descripiton area
        plot_description_layout = QtWdgt.QVBoxLayout()
        self.create_plot_description(plot_description_layout)
        user_input_result_layout.addLayout(plot_description_layout)
        user_input_result_layout.addStretch()
        
    # This function create the elements necessary for selecting the SANS data file
    def create_file_selection(self, file_selection_row):
        # Data file section
        file_selection_layout = QtWdgt.QHBoxLayout()

        # Button to open file dialog
        self.choose_file_dir_button = QtWdgt.QPushButton("Choose File")
        self.choose_file_dir_button.setFixedWidth(round(180*self.scale))
        self.choose_file_dir_button.setStyleSheet(
            f'''
            QPushButton{{padding: {round(6*self.scale)}px 0px {round(6*self.scale)}px 0px;
                         font-size: {round(18.5*self.scale)}px}}''')
        self.choose_file_dir_button.clicked.connect(self.choose_file)
        file_selection_layout.addWidget(self.choose_file_dir_button)
        
        # Tooltip
        file_select_description = (
            "The program reads ASCII data files. The first two required columns are the "
            "scattering vector Q and the scattering intensity I(Q), respectively.<br><br>"
            "An optional third column, dI(Q), can be included. If it cannot be acquired "
            "from the data file, percentage errors are to be used instead.<br><br>"
            "A wide range of delimiters is supported, but each file must use only one type of "
            "delimiter. Note: ',' <b>cannot</b> be used as a decimal point."
        )
        file_selection_tool_tip = self.make_tool_tip(file_select_description)   
        file_selection_layout.addWidget(file_selection_tool_tip)
        
        # Text box to display and edit the file path
        self.file_dir_label = QtWdgt.QLabel("Select data file...")
        self.file_dir_label.setStyleSheet(
            f'''
            QLabel{{font-size: {round(18*self.scale)}px;}}''')
        self.file_dir_label.setFixedWidth(self.input_result_width
                                          -round(180*self.scale))
        file_selection_layout.addWidget(self.file_dir_label)
        
        # Add element to the main program
        file_selection_row.addLayout(file_selection_layout)
        
    # This function create the elements necessary for changing the flat background and Q-max
    def create_bkgrd_Q_range_input(self, data_input_grid, row):
        # Background input elements
        bkgrd_label = QtWdgt.QLabel("Background (cm<sup>-1</sup>)")
        data_input_grid.addWidget(bkgrd_label, row, 0)
        self.bkgrd_input_box = QtWdgt.QLineEdit()
        self.bkgrd_input_box.setPlaceholderText("0.0")
        data_input_grid.addWidget(self.bkgrd_input_box, row, 
                                  self.input_field_col, 1, 2)
        bkgrd_description = (
            "<b>Background (cm<sup>-1</sup>):</b><br>"
            "This is the flat background value to be subtracted from the original SAS profile prior to the analysis.<br><br>"
            "The background value is determined such that the scattering profile forms a straight, continuous line in the large-Q region after subtraction."
            )
        bkgrd_tool_tip = self.make_tool_tip(bkgrd_description)
        data_input_grid.addWidget(bkgrd_tool_tip, row, self.input_tool_tip_col)
        
        # Q range input elements
        Q_range_label = QtWdgt.QLabel("Q range (\u212B<sup>-1</sup>)")
        data_input_grid.addWidget(Q_range_label, row+1, 0)
        self.Qmin_input_box = QtWdgt.QLineEdit()
        self.Qmin_input_box.setPlaceholderText("0")
        self.Qmax_input_box = QtWdgt.QLineEdit()
        self.Qmax_input_box.setPlaceholderText("\u221E")
        data_input_grid.addWidget(self.Qmin_input_box, row+1, 
                                  self.input_field_col)
        data_input_grid.addWidget(self.Qmax_input_box, row+1, 
                                  self.input_field_col+1)
        Qmax_description = (
            "<b>Q max (Å<sup>-1</sup>):</b> The maximum Q value to be used for the analysis."
            )
        Qmax_tool_tip = self.make_tool_tip(Qmax_description)
        data_input_grid.addWidget(Qmax_tool_tip, row+1, self.input_tool_tip_col)
        
        # Background input element
        dIQ_label = QtWdgt.QLabel("Measurement<br>error (cm<sup>-1</sup>)")
        data_input_grid.addWidget(dIQ_label, row+2, 0)
        
        ### Choose dIQ from data
        dIQ_data_layout = QtWdgt.QHBoxLayout()
        dIQ_data_button = QtWdgt.QRadioButton()
        dIQ_data_button.setStyleSheet(f'''
                                      QRadioButton{{padding-right: {round(-1.5*self.scale)}px;}}''')
        dIQ_data_label = QtWdgt.QLabel("From data")
        dIQ_data_layout.addStretch()
        dIQ_data_layout.addWidget(dIQ_data_button)
        dIQ_data_layout.addWidget(dIQ_data_label)
        dIQ_data_layout.addStretch()
        dIQ_data_layout.addStretch()
 
        ### Choose user-provided dIQ
        dIQ_user_layout = QtWdgt.QHBoxLayout()
        dIQ_user_button = QtWdgt.QRadioButton()
        dIQ_user_button.setStyleSheet(f'''
                                      QRadioButton{{padding-right: {round(-1.5*self.scale)}px;}}''')
        self.dIQ_percent_input_box = QtWdgt.QLineEdit()
        self.dIQ_percent_input_box.setPlaceholderText("2")
        self.dIQ_percent_input_box.setFixedWidth(round(30*self.scale))
        self.dIQ_percent_input_box.setAlignment(QtCore.Qt.AlignRight)
        dIQ_percent_label = QtWdgt.QLabel("% I(Q)")
        dIQ_user_layout.addStretch()
        dIQ_user_layout.addWidget(dIQ_user_button)
        dIQ_user_layout.addWidget(self.dIQ_percent_input_box)
        dIQ_user_layout.addWidget(dIQ_percent_label)
        dIQ_user_layout.addStretch()
        dIQ_user_layout.addStretch()
        
        ### Set button behaviour for selecting dIQ
        self.choose_dIQ = QtWdgt.QButtonGroup(self)
        self.choose_dIQ.setExclusive(True)  
        self.choose_dIQ.addButton(dIQ_data_button)
        self.choose_dIQ.addButton(dIQ_user_button)
        dIQ_data_button.setChecked(True)

        ### Combine elements and provide description
        data_input_grid.addLayout(dIQ_data_layout, row+2, 
                                  self.input_field_col)
        data_input_grid.addLayout(dIQ_user_layout, row+2, 
                                  self.input_field_col+1)
        dIQ_description = (
            "<b>dI(Q) (cm<sup>-1</sup>):</b> Measurement error of the scattering intensity. " + 
            "Can either be obtained from the input data file or specified as a percentage of I(Q)."
            )
        Qmax_tool_tip = self.make_tool_tip(dIQ_description)
        data_input_grid.addWidget(Qmax_tool_tip, row+2, self.input_tool_tip_col)
        
        # Confirm button
        self.confirm_bkgrd_Q_range_button = QtWdgt.QPushButton("Confirm\nValue")
        self.confirm_bkgrd_Q_range_button.setMinimumWidth(self.confirm_button_width)
        self.confirm_bkgrd_Q_range_button.setSizePolicy(QtWdgt.QSizePolicy.Policy.Expanding, 
                                                       QtWdgt.QSizePolicy.Policy.Expanding)
        self.confirm_bkgrd_Q_range_button.clicked.connect(self.set_bkgrd_Q_range)
        self.confirm_bkgrd_Q_range_button.setEnabled(False)
        data_input_grid.addWidget(self.confirm_bkgrd_Q_range_button, row, 
                                  self.input_confirm_col, 3, 1)
        
    # Create input fields for points-per-decade and smoothing factor λ,
    # and set up input monitoring via focus events and timers.
    def create_pts_per_dec_lambda_input(self, data_input_grid, row):
        def add_label_input_tooltip(label_text, placeholder, tooltip_html, target_attr, grid_row):
            """Helper to add label, QLineEdit, and tooltip to a grid row."""
            # Label
            label = QtWdgt.QLabel(label_text)
            data_input_grid.addWidget(label, grid_row, 0)
    
            # Input box
            input_box = QtWdgt.QLineEdit()
            input_box.setPlaceholderText(placeholder)
            data_input_grid.addWidget(input_box, grid_row, self.input_field_col, 1, 2)
    
            # Tooltip
            tool_tip = self.make_tool_tip(tooltip_html)
            data_input_grid.addWidget(tool_tip, grid_row, 4)
    
            # Store reference
            setattr(self, target_attr, input_box)
    
        # -- Add Points per Decade Input --
        pts_dec_label = "Number of points per\ndecade for result"
        pts_dec_tooltip = (
            "<b>Number of points per decade for result:</b><br>"
            "Refers to the number of r<sub>i</sub> values per decade of r, where r = 2.5/Q, in the fit results.<br><br>"
            "A higher number of r<sub>i</sub> values improves the fit to the original intensity profile. "
            "However, too many points may overfit the data and slow down computation.<br><br>"
            "Ensure the number of fit result points is less than half the number of input points."
        )
        add_label_input_tooltip(pts_dec_label, "10", pts_dec_tooltip, "pts_per_dec_input_box", row)
    
        # -- Add Lambda Input --
        lambda_label = "Smoothing factor λ"
        lambda_tooltip = (
            "<b>Smoothing factor (λ):</b><br>"
            "Controls the trade-off between fitting the data and smoothing the result—"
            "a higher λ leads to a smoother result.<br><br>"
            "The fitting process minimises Ξ = χ² + λ·ℜ, where χ² measures the data mismatch, and ℜ penalises roughness.<br><br>"
            "• λ = 0 disables smoothing (the fit follows the data closely).<br>"
            "• λ > 0 increases smoothness, possibly at the cost of a higher χ².<br><br>"
            "Try λ = 1 to start. Then increase or decrease it by factors of 10 to fine-tune."
            )        
        add_label_input_tooltip(lambda_label, "1", lambda_tooltip, "lambda_input_box", row + 1)
        
        # The following section help automatic detection of the value change in 
        # self.pts_per_dec_input_box and assign it to self.pts_per_dec
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(50)  # 1 second

        def make_focus_in_handler(widget):
            def focus_in(event):
                self.timer.start()
                QtWdgt.QLineEdit.focusInEvent(widget, event)
            return focus_in
        
        def make_focus_out_handler(widget):
            def focus_out(event):
                self.timer.stop()
                QtWdgt.QLineEdit.focusOutEvent(widget, event)
            return focus_out    
        
        def check_inputs():
            if not self.run_fit_button.isEnabled():
                # Points per decade
                try:
                    pts_text = self.pts_per_dec_input_box.text()
                    pts_val = 10 if pts_text == '' else int(pts_text)
                    if pts_val > 0 and pts_val != self.pts_per_dec:
                        if self.chosen_data_file_dir:
                            self.run_fit_button.setEnabled(True)
                except ValueError:
                    pass
                
                # Lambda
                try:
                    lam_text = self.lambda_input_box.text()
                    lam_val = 1 if lam_text == '' else float(lam_text)
                    if lam_val >= 0 and lam_val != self.lambda_:
                        if self.chosen_data_file_dir:
                            self.run_fit_button.setEnabled(True)
                except ValueError:
                    pass
            
        # Connect events
        self.pts_per_dec_input_box.focusInEvent = make_focus_in_handler(self.pts_per_dec_input_box)
        self.pts_per_dec_input_box.focusOutEvent = make_focus_out_handler(self.pts_per_dec_input_box)
        self.lambda_input_box.focusInEvent = make_focus_in_handler(self.lambda_input_box)
        self.lambda_input_box.focusOutEvent = make_focus_out_handler(self.lambda_input_box)
        self.timer.timeout.connect(check_inputs)
        
    # This function create elements in the program for handling the inputs required
    # for the fit of PDSP model
    def create_PDSP_fit_input(self, data_input_grid, row):
        # Input box for contrast between 2 phases
        contrast_label = QtWdgt.QLabel("Contrast between<br>2 phases "+
                                       "(10<sup>10</sup> cm<sup>-2</sup>)")
        data_input_grid.addWidget(contrast_label, row, 0)
        self.contrast_input_box = QtWdgt.QLineEdit()
        self.contrast_input_box.setPlaceholderText("3.0")
        data_input_grid.addWidget(self.contrast_input_box, row, 
                                  self.input_field_col, 1, 2)
        contrast_description = (
            "<b>Contrast between 2 phases (10<sup>10</sup> cm<sup>-2</sup>):</b><br>"
            "Refers to ρ<sub>1</sub><sup>*</sup> - ρ<sub>2</sub><sup>*</sup> in the PDSP model.<br>"
            "For porous systems without any filling medium (i.e., the filling medium is vacuum or air), the contrast is equal to the SLD of the solid."
            )
        contrast_tool_tip = self.make_tool_tip(contrast_description)
        data_input_grid.addWidget(contrast_tool_tip, row, self.input_tool_tip_col)
        
        # Input box for solid density
        density_label = QtWdgt.QLabel("Sample bulk density<br>(g/cm<sup>3</sup>)")
        data_input_grid.addWidget(density_label, row+1, 0)
        self.density_input_box = QtWdgt.QLineEdit()
        self.density_input_box.setPlaceholderText("1.0")
        data_input_grid.addWidget(self.density_input_box, row+1, 
                                  self.input_field_col, 1, 2)
        density_description = (
            "<b>Sample bulk density (g/cm<sup>3</sup>): </b><br>"
            "Determined as m<sub>sample</sub>/V<sub>sample</sub>, required for the volume-weighted pore distribution calculations."
            )
        density_tool_tip = self.make_tool_tip(density_description)
        data_input_grid.addWidget(density_tool_tip, row+1,
                                  self.input_tool_tip_col)
        
        # Input box for r value for SSA extrapolation
        r_SSA_extrapolate_label = QtWdgt.QLabel("Pore radius for SSA<br>extrapolation (nm)")
        data_input_grid.addWidget(r_SSA_extrapolate_label, row+2, 0)
        self.r_SSA_extrapolate_input_box = QtWdgt.QLineEdit()
        self.r_SSA_extrapolate_input_box.setPlaceholderText("0.2")
        data_input_grid.addWidget(self.r_SSA_extrapolate_input_box, row+2, 
                                  self.input_field_col, 1, 2)
        r_SSA_extrapolate_description = (
            "<b>Pore radius for SSA extrapolation (nm):</b><br>"
            "Is the pore radius to which the SSA value is extrapolated "
            "in the SSA(R) vs. R graph."
            )
        r_SSA_extrapolate_tool_tip = self.make_tool_tip(r_SSA_extrapolate_description)
        data_input_grid.addWidget(r_SSA_extrapolate_tool_tip, row+2,
                                  self.input_tool_tip_col)

        # Input box for number of points for SSA extrapolation
        num_pts_SSA_extrapolate_label = QtWdgt.QLabel("Number of points\nfor SSA extrapolation")
        data_input_grid.addWidget(num_pts_SSA_extrapolate_label, row+3, 0)
        self.num_pts_SSA_extrapolate_input_box = QtWdgt.QLineEdit()
        self.num_pts_SSA_extrapolate_input_box.setPlaceholderText("7")
        data_input_grid.addWidget(self.num_pts_SSA_extrapolate_input_box, row+3,
                                  self.input_field_col, 1, 2)
        num_pts_SSA_extrapolate_description = (
            "<b>Number of points for SSA extrapolation:</b><br>"
            "Refers to the number of data points used for the SSA extrapolation "
            "to the  specified pore radius in the SSA(R) vs. R graph."
            )
        num_pts_SSA_extrapolate_tool_tip = self.make_tool_tip(num_pts_SSA_extrapolate_description)
        data_input_grid.addWidget(num_pts_SSA_extrapolate_tool_tip, row+3,
                                  self.input_tool_tip_col)
        
        # Input major phase
        major_phase_label = QtWdgt.QLabel("Major phase")
        data_input_grid.addWidget(major_phase_label, row+4, 0)
        major_phase_description = (
            "<b>Major phase:</b><br>"
            "Since scattering intensity is symmetric for both phases, the phase "
            "occupying most of the system's volume must be selected to ensure the "
            "correct porosity value."
            )        
        major_phase_tool_tip = self.make_tool_tip(major_phase_description)
        data_input_grid.addWidget(major_phase_tool_tip, row+4, 
                                  self.input_tool_tip_col)
        
        ### Choose solid
        solid = QtWdgt.QHBoxLayout()
        solid_button = QtWdgt.QRadioButton('Solid')
        solid.addStretch()
        solid.addWidget(solid_button)
        solid.addStretch()
        data_input_grid.addLayout(solid, row+4, self.input_field_col)

        ### Choose void
        void = QtWdgt.QHBoxLayout()
        void_button = QtWdgt.QRadioButton('Void')
        void.addStretch()
        void.addWidget(void_button)
        void.addStretch()
        data_input_grid.addLayout(void, row+4, self.input_field_col+1)
        
        ### Set button behaviour for selecting dIQ
        self.choose_major_phase = QtWdgt.QButtonGroup(self)
        self.choose_major_phase.setExclusive(True)  
        self.choose_major_phase.addButton(solid_button)
        self.choose_major_phase.addButton(void_button)
        solid_button.setChecked(True)

        # Confirm button
        self.recalc_PDSP_input_button = QtWdgt.QPushButton("Recalculate\nFit Result")
        self.recalc_PDSP_input_button.setStyleSheet(f'''QPushButton{{font-weight: bold; font-size: {round(19*self.scale)}px;}}''')
        self.recalc_PDSP_input_button.setMinimumWidth(self.confirm_button_width)
        self.recalc_PDSP_input_button.setSizePolicy(QtWdgt.QSizePolicy.Policy.Expanding, 
                                                    QtWdgt.QSizePolicy.Policy.Expanding)
        self.recalc_PDSP_input_button.clicked.connect(self.recalc_PDSP_result)
        self.recalc_PDSP_input_button.setEnabled(False)
        data_input_grid.addWidget(self.recalc_PDSP_input_button, row, 
                                  self.input_confirm_col, 5, 1)
        
    # This function format the 'Run PDSP Fit' button
    def format_PDSP_button(self):
        self.run_fit_button.setStyleSheet(                                          
            f'''
            /* Enabled Button - Matching SP_MessageBoxQuestion */
            QPushButton {{
                font-size: {round(26*self.scale)}px;
                font: bold italic;
                color: white;
                padding: {round(12*self.scale)}px 0px {round(12*self.scale)}px 0px;
                background-color: #1565C0;  /* Deep Blue */
                border-radius: {round(8*self.scale)}px;
                }}
            /* Disabled Button - Reduced Contrast */
            QPushButton:disabled {{
                background-color: #64B5F6;  /* Lighter Blue */
                color: #90CAF9;  /* Even Lighter Blue for Text */
                }}
            '''
            )
        self.run_fit_button.clicked.connect(self.run_fit_func)
        self.run_fit_button.setEnabled(False)
        
    # This function create elements in the program for displaying the result
    # after the fit of PDSP model
    def create_result(self, result_grid):
        
        def add_result_entry(result_grid, label_text, row):
            label = QtWdgt.QLabel(label_text)
            result_grid.addWidget(label, row, 0)
            result_box = QtWdgt.QLineEdit()
            result_box.setReadOnly(True)
            result_box.setPlaceholderText('-')
            result_grid.addWidget(result_box, row, 2)
            
            return label, result_box  # Return both for future reference   
        
        _, self.porosity_box = add_result_entry(result_grid, 'Porosity', 0)
        _, self.pore_volume_box = \
            add_result_entry(result_grid, 
                             'Average pore<br>volume (cm<sup>3</sup>)', 1)
        _, self.pore_concentration_box = \
            add_result_entry(result_grid, 
                              'Pore concentration<br>(cm<sup>-3</sup>)', 2)        
        self.SSA_label, self.SSA_box = \
            add_result_entry(result_grid, 
                             ('SSA extrapolated<br>to {:.2f} nm'
                              .format(self.r_SSA_extrapolate) +
                              ' (cm<sup>2</sup>/cm<sup>3</sup>)'), 3)
        # Save result button
        self.save_result_button = QtWdgt.QPushButton('Save PDSP Result', self)
        self.save_result_button.setStyleSheet(
            f'''
            QPushButton{{font-size: {round(20*self.scale)}px; 
                         font: bold italic; 
                         padding: {round(7*self.scale)}px 0px {round(7*self.scale)}px 0px}}''')
        self.save_result_button.clicked.connect(self.save_result_func)
        self.save_result_button.setEnabled(False)
        result_grid.addWidget(self.save_result_button, 4, 0, 1, 3)
        
    # Add plot description for user
    def create_plot_description(self, plot_description_layout):
        description = QtWdgt.QLabel(
            f'''
            <div style='font-size: {round(self.scale*16)}px; 
                font-style: italic; margin: {round(10*self.scale)}px;'>
            <p style='margin-bottom: {5*self.scale}px;'>
            <b>dV/dr: </b>
            Differential pore volume per unit weight.</p>

            <p style='margin-bottom: {5*self.scale}px;'>
            <b>f(r): </b>
            Probability density function of the pore sizes (number-weighted pore size distribution).</p>

            <p>
            <b>SSA(R): </b>
            Specific surface area for probe with radius R, 
            (defined as the sum of SSA's of all pores with radii larger than R, 
                 divided by the sample volume).</p>
             </div>
             ''')
        description.setWordWrap(True)
        plot_description_layout.addWidget(description)

    # Creating 4 plotting windows for visualising PDSP fit result
    def create_figure_windows(self):
        # Plotting window for SAS data and background-subtracted SAS data
        plot_group_box_SAS, self.figure_SAS, self.canvas_SAS =\
                                        self.create_plot_layout('SAS Data')
        pf.set_SAS_plot(self.figure_SAS, self.scale)
        self.main_layout.addWidget(plot_group_box_SAS, 0, 1)
        
        # Plotting window for background-subtracted data and fitted data
        plot_group_box_SAS_fitted, self.figure_SAS_fitted, self.canvas_SAS_fitted  =\
                        self.create_plot_layout('SAS Data vs. Fitted Result')
        pf.set_SAS_plot(self.figure_SAS_fitted, self.scale)
        self.main_layout.addWidget(plot_group_box_SAS_fitted, 0, 2)
        
        # Plotting window for dV/dr vs. r plot, resulting from the fit of the PDSP model
        plot_group_box_dVdr, self.figure_dVdr, self.canvas_dVdr =\
                                        self.create_plot_layout('dV/dr Plot')
        pf.set_dVdr_plot(self.figure_dVdr, self.scale)
        self.main_layout.addWidget(plot_group_box_dVdr, 1, 1)

        # Plotting window for f(r) vs. r and SSA(R) vs R plot, 
        # resulting from the fit of the PDSP model
        plot_group_box_fr_SSA, self.figure_fr_SSA, self.canvas_fr_SSA =\
                            self.create_plot_layout('f(r) vs. r  ||  SSA(R) vs. R')
        pf.set_fr_SSA_plot(self.figure_fr_SSA, self.scale)
        self.main_layout.addWidget(plot_group_box_fr_SSA, 1, 2)
        
    # This function used to create the individual plotting window
    def create_plot_layout(self, plot_name):
        group_box = QtWdgt.QGroupBox(plot_name)
        plot_layout = QtWdgt.QVBoxLayout()
        plot_layout.setSpacing(round(6*self.scale))
        group_box.setStyleSheet(
            f'''
            QGroupBox{{font: bold {round(19*self.scale)}px Arial;  
                       border: {round(1*self.scale)}px solid; 
                       margin-top: {round(15*self.scale)}px; 
                       padding-top: {round(10*self.scale)}px;}}''' + 
            f'''
            QGroupBox:title{{subcontrol-origin: margin;
                             padding: 0px {round(5*self.scale)}px 0px {round(5*self.scale)}px; 
                             subcontrol-position: top center;}}''')
        fig = mpl_figure.Figure() # assign new figure
        canvas = mpl_backend.FigureCanvas(fig) # assign new plotting canvas
        canvas.setMinimumSize(round(600*self.scale), round(525*self.scale))
        toolbar = mpl_backend.NavigationToolbar2QT(canvas, self) # add figure tool bar
        toolbar.setStyleSheet(
            f'''
            font-size:{round(14*self.scale)}px;''')
        toolbar.setStyleSheet(
            f'''
            QToolButton{{ 
                width:{round(26*self.scale)}px;
                height:{round(26*self.scale)}px;
                }}''')
        plot_layout.addWidget(toolbar)
        plot_layout.addWidget(canvas)
        group_box.setLayout(plot_layout)
        return group_box, fig, canvas

    # Function create for selecting the SAS data file, activated when the button
    # 'Choose File' is clicked
    def choose_file(self):
        # Open file dialog
        file_dir, _ = QtWdgt.QFileDialog.getOpenFileName(self, 
                                                         "Select SANS Data File", 
                                                         self.chosen_data_folder_dir or "", 
                                                         "Text Files (*.txt *.dat *.csv *.ABS)")
        if file_dir:  # If a file is selected
            # Read and obtain data from the selected file, then assign the
            # data to the corresponding predefined variables and 
            try:
                self.QQ_origin, self.IQ_origin, self.dIQ_data = \
                                                    bf.read_SANS_data(file_dir)
            except ValueError as e:
                self.show_error_message(str(e))
                return
            
            if isinstance(self.dIQ_data, int) or \
                    (isinstance(self.dIQ_data, (list, np.ndarray)) and (len(self.dIQ_data) != len(self.QQ_origin) or 
                                                                        np.sum(self.dIQ_data/self.IQ_origin) < 1e-5)):
                self.choose_dIQ.buttons()[0].setEnabled(False)
                self.choose_dIQ.buttons()[0].setChecked(False)
                self.choose_dIQ.buttons()[1].setChecked(True)
            else:
                self.choose_dIQ.buttons()[0].setEnabled(True)

            if self.choose_dIQ.buttons()[1].isChecked():
                try:
                    self.set_parameter("IQ_percent_dIQ", 'Percentage Error',self.dIQ_percent_input_box, 
                                       default_val = 2, min_val = 0, 
                                       error_msg = ("Invalid dI(Q) percentage value. Must be a number > 0. " +
                                                    "Background percentage set to 2%"))
                except ValueError as e:
                    self.show_error_message(str(e))
                self.dIQ_user = self.IQ_origin * self.IQ_percent_dIQ/100

            # Save chosen file location and displaying file info
            self.chosen_data_file_dir = file_dir
            self.chosen_data_folder_dir = '/'.join(self.chosen_data_file_dir.split('/')[:-1])
            self.file_dir_label.setText(self.chosen_data_file_dir.split('/')[-1]) # Display selected file name
            print('Folder of chosen data file: ' + self.chosen_data_folder_dir)
            print('')
            
            # Clear plots, result prior to plotting new data set
            bf.clear_plot(self.figure_SAS, self.canvas_SAS)
            bf.clear_plot(self.figure_SAS_fitted, self.canvas_SAS_fitted)
            self.clear_result()
                        
            # Plot the input SAS data
            self.dIQ_origin = self.dIQ_data.copy() \
                if self.choose_dIQ.buttons()[0].isChecked() \
                    else self.dIQ_user.copy()
                
            bf.plot_SANS_data(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                              self.figure_SAS, self.canvas_SAS)
            self.QQ_trim = self.QQ_origin.copy()
            self.IQ_trim = self.IQ_origin.copy()
            self.dIQ_trim = self.dIQ_origin.copy()
            self.QQ_trim, self.IQ_trim, self.dIQ_trim = \
                bf.subtract_background(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                                       self.bkgrd, self.Qmin, self.Qmax)
            bf.plot_SANS_subtract(self.QQ_trim, self.IQ_trim, 
                                  self.QQ_origin, self.bkgrd,
                                  self.figure_SAS, self.canvas_SAS)
            bf.plot_SANS_fit(self.QQ_trim, self.IQ_trim, 
                             self.figure_SAS_fitted, self.canvas_SAS_fitted, 
                             which = 'input')
            
            # Enable/disable action buttons to prevent accidental inputs
            self.confirm_bkgrd_Q_range_button.setEnabled(True)
            self.run_fit_button.setEnabled(True)       
            self.recalc_PDSP_input_button.setEnabled(False)
            self.save_result_button.setEnabled(False)


    # Function to assign value to fitting variables with validation.
    def set_parameter(self, attr, attr_descrpition, input_box, 
                      default_val, min_val, error_msg):
        # List of attributes that allow values to be greater than 
        # or equal to (>=) the minimum value.
        attr_larger_equal = ['bkgrd', 'Qmin', 'pts_per_dec',
                             'num_pts_SSA_extrapolate',]
        text = input_box.text()
        if text:
            try:
                # Convert the input text to the appropriate numeric type
                # based on the parameter.
                value = (int(text) if (attr == 'num_pts_SSA_extrapolate' or
                                       attr == 'pts_per_dec')
                         else  float(text)*1e10 if attr == 'contrast'
                         else float(text))
                
                # Check if the value meets the minimum requirement
                # and assign it to the attribute.
                if (value >= min_val 
                    if attr in attr_larger_equal else value > min_val):       
                    setattr(self, attr, value)
                    print(f"{attr_descrpition} set to: {value}" + 
                          ("%" if 'percent' in attr else ""))
                else:
                    raise ValueError(error_msg)
            except ValueError:
                setattr(self, attr, default_val)
                input_box.setText('')
                print(f"{attr_descrpition} set to: {default_val}" + 
                      ("%" if 'percent' in attr else ""))
                raise ValueError(error_msg)
        else:
            setattr(self, attr, default_val)
            print(f"{attr_descrpition} set to: {default_val}" + 
                  ("%" if 'percent' in attr else ""))
            
    # Function create for confirming the background being subtracted and the 
    # limiting Q range, activated when the button 'Confirm Value' is clicked
    def set_bkgrd_Q_range(self, clear_result = False, propagate_error = False):
        prev_values = {"bkgrd": self.bkgrd, "Qmin": self.Qmin, 
                       "Qmax": self.Qmax, "dIQ_origin": self.dIQ_origin}
        button_dIQ_data, button_dIQ_user = self.choose_dIQ.buttons()

        try:
            # Process Background, Qmin, and Qmax Inputs
            self.set_parameter("bkgrd", 'Background value',self.bkgrd_input_box, 
                               default_val = 0, min_val = 0, 
                               error_msg = "Invalid background value. Must be a number >= 0.")
            self.set_parameter("Qmin", 'Minimum Q value', self.Qmin_input_box,
                               default_val = 0, min_val=0, 
                               error_msg = "Invalid Q min value. Must be a number >= 0.")
            self.set_parameter("Qmax", 'Maximum Q value', self.Qmax_input_box, 
                               default_val = np.inf, min_val=0, 
                               error_msg = "Invalid Q max value. Must be a number > 0.")
            if button_dIQ_data.isChecked():
                self.dIQ_origin = self.dIQ_data.copy()
                print("dI(Q) error set to: From data")
            elif button_dIQ_user.isChecked():
                self.set_parameter("IQ_percent_dIQ", 'dI(Q) percentage error', self.dIQ_percent_input_box, 
                                   default_val = 2, min_val = 0, 
                                   error_msg = ("Invalid dI(Q) percentage value. Must be a number > 0. " +
                                                "dI(Q) percentage set to 2%"))
                self.dIQ_user = self.IQ_origin * self.IQ_percent_dIQ/100
                self.dIQ_origin = self.dIQ_user.copy()
                
        except ValueError as e:
            if not propagate_error:
                self.show_error_message(str(e))
                return
            else:
                raise ValueError(str(e))
        
        print('')
        value_changed = any(not np.array_equal(prev_values[attr], getattr(self, attr))
                            if isinstance(prev_values[attr], np.ndarray) 
                            else prev_values[attr] != getattr(self, attr)
                            for attr in prev_values)

        # Replot subtracted data if any values changed
        if value_changed or clear_result:            
            if not np.array_equal(prev_values['dIQ_origin'], self.dIQ_origin):
                self.QQ_trim, self.IQ_trim, self.dIQ_trim = \
                    bf.subtract_background(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                                           self.bkgrd, self.Qmin, self.Qmax)
                bf.plot_SANS_data(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                                  self.figure_SAS, self.canvas_SAS)

            self.QQ_trim, self.IQ_trim, self.dIQ_trim = \
                bf.subtract_background(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                                       self.bkgrd, self.Qmin, self.Qmax)

            bf.plot_SANS_subtract(self.QQ_trim, self.IQ_trim, self.QQ_origin, self.bkgrd, self.figure_SAS, self.canvas_SAS)
            bf.plot_SANS_fit(self.QQ_trim, self.IQ_trim, self.figure_SAS_fitted, self.canvas_SAS_fitted, which='input')
            self.clear_result()
            self.recalc_PDSP_input_button.setEnabled(False)
            self.run_fit_button.setEnabled(True)
            self.save_result_button.setEnabled(False)
            
    # Set PDSP fit parameters: contrast, bulk density, pore radius for SSA extrapolate
    # number of points for SSA extrapolation, and major phase
    def set_PDSP_fit_inputs(self, propagate_error = False):
        try:
            self.set_parameter('contrast', 'Contrast between 2 phases', 
                               self.contrast_input_box, default_val = 3e10, min_val = 0, 
                               error_msg = 'Invalid contrast value. Must be a number > 0')
            self.set_parameter('density', 'Sample bulk density', 
                               self.density_input_box, default_val = 1.0, min_val = 0, 
                               error_msg = 'Invalid density value. Must be a number > 0')
            self.set_parameter('num_pts_SSA_extrapolate', 
                               'Number of points for for SSA extrapolation', 
                               self.num_pts_SSA_extrapolate_input_box, 
                               default_val = 7, min_val = 3, 
                               error_msg = ('Invalid number of points for SSA extrapolation. '+
                                            'Must be an integer >= 3.'))
            prev_val = self.r_SSA_extrapolate
            self.set_parameter('r_SSA_extrapolate', 'Pore radius for SSA extrapolation',
                               self.r_SSA_extrapolate_input_box, 
                               default_val = 0.2, min_val = 0,
                               error_msg = ('Invalid radius value for SSA extrapolation. ' +
                                            'Must be a number > 0'))
            if abs(prev_val - self.r_SSA_extrapolate) > 1e-8:
                self.result_grid.removeWidget(self.SSA_label)
                self.SSA_label = QtWdgt.QLabel('SSA extrapolated<br>to {:.2f} nm'
                                               .format(self.r_SSA_extrapolate) +
                                               ' (cm<sup>2</sup>/cm<sup>3</sup>)')
                self.SSA_label.setStyleSheet(f'''QLabel{{font-size: {round(19*self.scale)}px; }}''')
                self.result_grid.addWidget(self.SSA_label, 3, 0)
        except ValueError as e:
            if not propagate_error:
                self.show_error_message(str(e))
                return
            else:
                raise ValueError(str(e))
                
        self.major_phase = self.choose_major_phase.checkedButton().text().lower()
        print(f"Major phase of sample: {self.major_phase}")
        print('')
                
    # This function fits, returns, and plots the PDSP result based on the provided
    # inputs. Activated when the button 'Fit PDSP Model!' is clicked.
    def run_fit_func(self):
        try:
            # Confirming fitting inputs
            self.set_bkgrd_Q_range(clear_result=True, propagate_error = True)
            self.set_parameter('pts_per_dec', 'Number of points per decade for result', 
                               self.pts_per_dec_input_box, default_val = 10, min_val = 3, 
                               error_msg = ('Invalid number of points per decade for result. '
                                            'Must be an integer >= 3'))
            self.set_parameter('lambda_', 'Smoothing factor lambda', 
                               self.lambda_input_box, default_val = 1, min_val = 0, 
                               error_msg = ('Invalid lambda value. Must be > 0'))
            print('')
            self.set_PDSP_fit_inputs(propagate_error = True)
        except ValueError as e:
            self.show_error_message(str(e))
            return
            
        # Update program to show old results have been cleared
        QtWdgt.QApplication.processEvents()
        time.sleep(0.001)  
        
        # Run fit function
        print('Start PDSP fit')
        try:
            self.rr, self.IQ_fitted, self.IQ0_fitted, self.f_r, self.f_dash_r, self.SSA,\
                self.dV_dr, self.phi, self.Vpore_avg, self.phi_on_Vavg, self.SSA_extrapolate = \
                    bf.fit_PDSP_model(self.QQ_trim, self.IQ_trim, self.dIQ_trim,
                                      self.pts_per_dec, self.lambda_, 
                                      self.contrast, self.density, self.r_SSA_extrapolate, 
                                      self.num_pts_SSA_extrapolate, self.major_phase)
        except ValueError as e:
            self.show_error_message(str(e))
            return
        print('Done PDSP fit!\n')
        
        # Display result, enable/disable corresponding buttons to prevent accidental inputs.
        self.display_result()        
        self.recalc_PDSP_input_button.setEnabled(True)
        self.run_fit_button.setEnabled(False)
        self.save_result_button.setEnabled(True)
        
        # Return warning message if result are not real value
        if np.iscomplex(self.phi[0]):
            warning_message = (
                'Resulting <i>\u03D5</i>(1 - <i>\u03D5</i>) > 0.25, '
                'PDSP fit results are not real.<br>'
                'Check the absolute scattering intensity or contrast between two phases.')
            self.show_error_message(warning_message, error_type = 'Fit Error')
            
    # Recalculate PDSP results that doesn't require a full refit, activated
    # when the button 'Recalculate Fit Result' is pressed.
    def recalc_PDSP_result(self):
        # Confirm fitting parameters
        try:
            self.set_PDSP_fit_inputs(propagate_error = True)
        except ValueError as e:
            self.show_error_message(str(e))
            return

        # Recalculate fit result and display the result
        _, self.SSA, self.dV_dr, self.phi, self.Vpore_avg,\
            self.phi_on_Vavg, self.SSA_extrapolate =\
                bf.calc_PDSP_result(self.rr*10, self.f_r, self.f_dash_r, self.IQ0_fitted, 
                                    self.contrast, self.density, self.r_SSA_extrapolate,
                                    self.num_pts_SSA_extrapolate, self.major_phase)
        self.display_result()
        
    # Function used to present the result, including filling out the text boxes
    # in the result area and plotting the result.
    def display_result(self):
        # Filling out the text boxes in the result area
        self.porosity_box.setText('{:.3f} ± {:.3f}'.format(self.phi[0], self.phi[1]*self.phi[0]))
        self.pore_concentration_box.setText('{:.3e} ± {:.3e}'.format(self.phi_on_Vavg[0], 
                                                                     self.phi_on_Vavg[1]*self.phi_on_Vavg[0]))
        self.pore_volume_box.setText('{:.3e} ± {:.3e}'.format(self.Vpore_avg[0],
                                                              self.Vpore_avg[1]*self.Vpore_avg[0]))
        self.SSA_box.setText('{:.3e} ± {:.3e}'.format(self.SSA_extrapolate[0],
                                                      self.SSA_extrapolate[1]*self.SSA_extrapolate[0]))
        
        # Plotting result
        bf.plot_SANS_fit(self.QQ_trim, self.IQ_fitted, self.figure_SAS_fitted, 
                         self.canvas_SAS_fitted, which = 'result')      
        bf.plot_dVdr(self.rr, self.dV_dr, self.figure_dVdr, self.canvas_dVdr)   
        bf.plot_fr_SSA(self.rr, self.f_r, self.SSA, self.num_pts_SSA_extrapolate, 
                       self.r_SSA_extrapolate, self.figure_fr_SSA, self.canvas_fr_SSA)
                
    # Function used to clear the result, including clearing out the text boxes
    # in the result area and removing the plot of the PDSP fit result.
    def clear_result(self):
        self.porosity_box.setText('')
        self.pore_concentration_box.setText('')
        self.pore_volume_box.setText('')
        self.SSA_box.setText('')
        bf.clear_plot(self.figure_dVdr, self.canvas_dVdr)
        bf.clear_plot(self.figure_fr_SSA, self.canvas_fr_SSA)
    
    # Function used to save the PDSP fit result into a text file, activated
    # when the button 'Save PDSP Result' is pressed.
    def save_result_func(self):
        # create result table for r vs f(r), SSA(R), and dV/dr 
        data_table = np.column_stack((self.rr, self.f_r[:,0], self.SSA[:,0], self.dV_dr[:,0], self.IQ0_fitted[:,1]))
        
        # initialise save location
        file_name_save = (self.chosen_data_file_dir.split('/')[-1].replace('.txt','').replace('.ABS','')
                                  .replace('.dat','').replace('.csv','') 
                                  + " PDSP Result.txt")
        if self.chosen_save_folder_dir == '':
            self.chosen_save_folder_dir = self.chosen_data_folder_dir[:]
        save_file_dir_default = self.chosen_save_folder_dir + '/' + file_name_save
        
        # Open file dialog to save file
        options = QtWdgt.QFileDialog.Options()
        save_file_dir, _ = QtWdgt.QFileDialog.getSaveFileName(self, "Save File", 
                                                              save_file_dir_default, 
                                                              "Text Files (*.txt);;All Files (*)", options=options)
        # Remember previously used save folder
        self.chosen_save_folder_dir = '/'.join(save_file_dir.split('/')[:-1])
        
        if save_file_dir:
            # Create the file and write result
            with open(save_file_dir, 'w') as file:
                # Write file header, including file name, background value,
                # Q-max, contrast, solid density, porosity, average pore volume,
                # pore concentration, extrapolated SSA
                file.write('PDSP Fit Result for ' + self.chosen_data_file_dir.split('/')[-1])
                file.write('\n\n')
                file.write('Background value (cm-1): {:.3e}'.format(self.bkgrd))
                file.write('\n')
                file.write('Selected Q range (A-1): [{:.3e}, {:.3e}]'.format(self.Qmin, self.Qmax))
                file.write('\n')
                file.write('Smoothing factor Lambda: {:.1e}'.format(self.lambda_))
                file.write('\n')
                file.write('Contrast between 2 phases (cm-2): {:.3e}'.format(self.contrast))
                file.write('\n')
                file.write('Density of Solid (g/cm3): {:.3f}'.format(self.density))
                file.write('\n\n')            
                file.write('Porosity: {:.5e} ± {:.5e}'.format(self.phi[0], self.phi[1]*self.phi[0]))
                file.write('\n')
                file.write('Average Pore Volume (cm3): {:.3e} ± {:.3e}'.format(self.Vpore_avg[0], 
                                                                               self.Vpore_avg[1]*self.Vpore_avg[0]))
                file.write('\n')
                file.write('Pore Concentration (cm-3): {:.3e} ± {:.3e}'.format(self.phi_on_Vavg[0], 
                                                                               self.phi_on_Vavg[1]*self.phi_on_Vavg[0]))
                file.write('\n')
                file.write('SSA interpolated to r = {:.2f} nm (cm2/cm3): '
                           .format(self.r_SSA_extrapolate) +
                           '{:.3e} ± {:.3e}'.format(self.SSA_extrapolate[0],
                                                    self.SSA_extrapolate[1]*self.SSA_extrapolate[0]))
                file.write('\n\n')
                
                # Write result table for r vs f(r), SSA(R), and dV/dr
                file.write('Pore size distribution table\n')
                file.write('\t\t'.join(['r\t', 'f(r)', 'SSA\t', 'dV/dr', '% error']) + '\n')
                [file.write('\t'.join(val if isinstance(val, str)
                                      else '{:.5e}'.format(val) 
                                      if val > 0 
                                      else '{:.4e}'.format(val) 
                                      for val in line) + '\n')
                 for line in data_table]
            print(f"File saved as {file_name_save} in {self.chosen_save_folder_dir}\n")

    # Function used for creating tool tips for inputs and results
    def make_tool_tip(self, tool_tip_message):
        # Create tool tip icon
        tool_tip = QtWdgt.QLabel()
        tool_tip.setAlignment(QtCore.Qt.AlignCenter)  # Center align icon
        icon = (QtWdgt.QApplication.style()             # get PyQt builtin icon
                .standardIcon(QtWdgt.QApplication
                              .style().SP_MessageBoxQuestion))
        tool_tip.setPixmap(icon.pixmap(int(25*self.scale), 
                                       int(25*self.scale)))  # Set icon size
        tool_tip.setToolTip(tool_tip_message)
        return tool_tip
    
    # Show input of error without crashing PyInstaller
    def show_error_message(self, message, error_type = 'Input Error'):
        msg_box = QtWdgt.QMessageBox()
        
        # Determine the appropriate icon
        if error_type == 'Input Error':
            icon = QtWdgt.QMessageBox.Critical
            window_icon = QtWdgt.QApplication.style().standardIcon(QtWdgt.QStyle.StandardPixmap.SP_MessageBoxCritical)
        else:
            icon = QtWdgt.QMessageBox.Warning
            window_icon = QtWdgt.QApplication.style().standardIcon(QtWdgt.QStyle.StandardPixmap.SP_MessageBoxWarning)
        
        # Set the message box
        msg_box.setWindowIcon(QtGui.QIcon(window_icon))
        msg_box.setIcon(icon)
        msg_box.setWindowTitle(error_type)
        msg_box.setText(message)
        msg_box.setStandardButtons(QtWdgt.QMessageBox.Ok)
    
        msg_box.exec_()
            
if __name__ == "__main__":
    app = QtWdgt.QApplication(sys.argv)
    main_window = PRINSAS_App()
    main_window.show()
    sys.exit(app.exec_())