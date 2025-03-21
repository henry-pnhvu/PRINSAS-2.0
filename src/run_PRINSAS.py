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
        print('UI Scale factor = {:.2f}'.format(self.scale))
        
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
        self.section_spacing = round(20*self.scale)
        self.input_result_width = round(700*self.scale)
        self.confirm_button_width = round(180*self.scale)
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
        
        # SANS data after background subtraction
        self.QQ_subtract = []
        self.IQ_subtract = []
        self.dIQ_subtract = []
        
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
        # Set universal font to Arial with the size of 20px
        self.setStyleSheet(f'''QWidget{{font: Arial; font-size: {round(self.scale*20)}px;}}''')
        
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
        self.set_style_sheet_group_box(input_group_box, font_size=24, 
                                       groupbox_margin_top=15, 
                                       groupbox_padding_top=20, 
                                       title_padding_left=30)
        data_input_grid = QtWdgt.QGridLayout()
        data_input_grid.setRowMinimumHeight(3, round(30*self.scale))
        data_input_grid.setColumnMinimumWidth(1,round(17*self.scale))
        data_input_grid.setColumnMinimumWidth(2,round(120*self.scale))
        data_input_grid.setColumnMinimumWidth(3,round(40*self.scale))
        data_input_grid.setVerticalSpacing(round(15*self.scale))
        input_group_box.setLayout(data_input_grid)
        user_input_result_layout.addWidget(input_group_box)
        user_input_result_layout.insertSpacing(1, self.section_spacing)
        
        # Add background input row to user input area
        self.create_bkgrd_Q_range_input(data_input_grid, row = 0)
        
        # Add number of points per decade input row to user input area
        self.create_pts_per_dec_input(data_input_grid, row = 2)
        
        # Add PDSP input row to user input area
        self.create_PDSP_fit_input(data_input_grid, row = 4)
        
        # Add run fit button to user input area
        self.run_fit_button = QtWdgt.QPushButton("Fit PDSP Model!")
        self.format_PDSP_button()
        user_input_result_layout.addWidget(self.run_fit_button)
        user_input_result_layout.insertSpacing(3, self.section_spacing)

        # Add result display area
        result_group_box = QtWdgt.QGroupBox("PDSP Fit Results")
        self.set_style_sheet_group_box(result_group_box, font_size=24, 
                                       groupbox_margin_top=15, 
                                       groupbox_padding_top=20, 
                                       title_padding_left=30)
        user_input_result_layout.addWidget(result_group_box)
        user_input_result_layout.insertSpacing(5, self.section_spacing)
        self.result_grid = QtWdgt.QGridLayout()
        self.result_grid.setColumnMinimumWidth(1,round(5*self.scale))
        self.result_grid.setVerticalSpacing(round(15*self.scale))
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
            QPushButton{{padding: {round(6*self.scale)}px 0px {round(6*self.scale)}px 0px;}}''')
        self.choose_file_dir_button.clicked.connect(self.choose_file)
        file_selection_layout.addWidget(self.choose_file_dir_button)

        # Text box to display and edit the file path
        self.file_dir_label = QtWdgt.QLabel("Select data file...")
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
            "This is the flat background value to be subtracted from the original SAS profile prior to the analysis.<br>"
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
        
        # Confirm button
        self.confirm_bkgrd_Q_range_button = QtWdgt.QPushButton("Confirm\nValue")
        self.confirm_bkgrd_Q_range_button.setMinimumWidth(self.confirm_button_width)
        self.confirm_bkgrd_Q_range_button.setSizePolicy(QtWdgt.QSizePolicy.Policy.Expanding, 
                                                       QtWdgt.QSizePolicy.Policy.Expanding)
        self.confirm_bkgrd_Q_range_button.clicked.connect(self.set_bkgrd_Q_range)
        self.confirm_bkgrd_Q_range_button.setEnabled(False)
        data_input_grid.addWidget(self.confirm_bkgrd_Q_range_button, row, 
                                  self.input_confirm_col, 2, 1)
        
    # This function create the elements necessary for changing the number of points 
    # per decade of the fit result
    def create_pts_per_dec_input(self, data_input_grid, row):
        # Label
        pts_per_dec_label = QtWdgt.QLabel("Number of points per\n" +
                                          "decade for result")
        data_input_grid.addWidget(pts_per_dec_label, row, 0)
        
        # Text box for inputing Pts per dec value
        self.pts_per_dec_input_box = QtWdgt.QLineEdit()
        self.pts_per_dec_input_box.setPlaceholderText("10")
        data_input_grid.addWidget(self.pts_per_dec_input_box, row, 
                                  self.input_field_col, 1, 2)
        # Tool tip element
        description = (
            "<b>Number of points per decade for result:</b><br>"
            "Refers to the number of r<sub>i</sub> values per decade of r, where r = 2.5/Q, in the fit results.<br>"
            "A higher number of r<sub>i</sub> values improves the fit to the original intensity profile. However, too many r<sub>i</sub> values lead to overfitting and increase computation time.<br>"
            "Should be chosen such that the number of points in the fit result remains fewer than half of the input data points."
            )       
        tool_tip = self.make_tool_tip(description)
        data_input_grid.addWidget(tool_tip, row, 4)
        
        # The following section help automatic detection of the value change in 
        # self.pts_per_dec_input_box and assign it to self.pts_per_dec
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(50)  # 1 second

        # Focus in event handler
        def on_focus_in(event):
            """Handle focus in event."""
            print("pts_per_dec_input_box focused")
            self.timer.start()  # Start the timer when focused
            super(QtWdgt.QLineEdit, self.pts_per_dec_input_box).focusInEvent(event)  # Call the base class method
        
        # Focus out event handler
        def on_focus_out(event):
            """Handle focus out event."""
            print("pts_per_dec_input_box lost focus")
            self.timer.stop()  # Stop the timer when focus is lost
            super(QtWdgt.QLineEdit, self.pts_per_dec_input_box).focusOutEvent(event)  # Call the base class method
        
        # Check input value on timer
        def check_pts_per_dec_input():
            """Check and print the current value of pts_per_dec_input_box."""
            try:
                if self.pts_per_dec_input_box.text() == '':
                    current_value = 10  # Default value
                else:
                    current_value = int(self.pts_per_dec_input_box.text())
                if current_value != self.pts_per_dec:
                    if current_value <= 0:
                        print("Invalid number of points per decade value. Please enter a positive integer.")
                    else:
                        self.pts_per_dec = current_value
                        print(f"Number of points per decade set to: {self.pts_per_dec}")
                        if self.chosen_data_file_dir != '':
                            self.run_fit_button.setEnabled(True)
            except ValueError:
                print("Invalid number of points per decade value. Please enter a positive integer.")
        
        self.pts_per_dec_input_box.focusInEvent = on_focus_in
        self.pts_per_dec_input_box.focusOutEvent = on_focus_out
        self.timer.timeout.connect(check_pts_per_dec_input)
        
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
        # Flip switch
        selection_switch = QtWdgt.QHBoxLayout()
        solid_label = QtWdgt.QLabel("Solid")
        void_label = QtWdgt.QLabel("Void")
        phase_switch = LeftRightSwitch(self)
        selection_switch.addWidget(solid_label, alignment=QtCore.Qt.AlignCenter)
        selection_switch.addWidget(phase_switch, alignment=QtCore.Qt.AlignCenter)
        selection_switch.addWidget(void_label, alignment=QtCore.Qt.AlignCenter)
        data_input_grid.addLayout(selection_switch, row+4, 
                                  self.input_field_col, 1, 2)

        # Confirm button
        self.recalc_PDSP_input_button = QtWdgt.QPushButton("Recalculate\nFit Result")
        self.recalc_PDSP_input_button.setStyleSheet(f'''QPushButton{{font-size: {round(25*self.scale)}px;}}''')
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
            <div style='font-size: {round(self.scale*18)}px; 
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
                self.QQ_origin, self.IQ_origin, self.dIQ_origin = \
                                                    bf.read_SANS_data(file_dir)
            except ValueError as e:
                self.show_error_message(str(e))
                return
            
            # Save chosen file location and displaying file info
            self.chosen_data_file_dir = file_dir
            self.chosen_data_folder_dir = '/'.join(self.chosen_data_file_dir.split('/')[:-1])
            self.file_dir_label.setText(self.chosen_data_file_dir.split('/')[-1]) # Display selected file name
            print('Folder of chosen data file: ' + self.chosen_data_folder_dir)
            
            # Clear plots, result prior to plotting new data set
            bf.clear_plot(self.figure_SAS, self.canvas_SAS)
            bf.clear_plot(self.figure_SAS_fitted, self.canvas_SAS_fitted)
            self.clear_result()
                        
            # Plot the input SAS data
            bf.plot_SANS_data(self.QQ_origin, self.IQ_origin, self.dIQ_origin,
                              self.figure_SAS, self.canvas_SAS)
            self.QQ_subtract = self.QQ_origin.copy()
            self.QQ_trim = self.QQ_origin.copy()
            self.IQ_subtract = self.IQ_origin.copy()
            self.IQ_trim = self.IQ_origin.copy()
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
        attr_larger_equal = ['bkgrd', 'Qmin', 'num_pts_SSA_extrapolate']
        text = input_box.text()
        if text:
            try:
                # Convert the input text to the appropriate numeric type
                # based on the parameter.
                value = (int(text) if attr == 'num_pts_SSA_extrapolate' 
                         else  float(text)*1e10 if attr == 'contrast'
                         else float(text))
                
                # Check if the value meets the minimum requirement
                # and assign it to the attribute.
                if (value >= min_val 
                    if attr in attr_larger_equal else value > min_val):       
                    setattr(self, attr, value)
                    print(f"{attr_descrpition} set to: {value}")
                else:
                    raise ValueError(error_msg)
            except ValueError:
                raise ValueError(error_msg)
        else:
            setattr(self, attr, default_val)
            print(f"{attr_descrpition} set to: {default_val}")
            
    # Function create for confirming the background being subtracted and the 
    # limiting Q range, activated when the button 'Confirm Value' is clicked
    def set_bkgrd_Q_range(self, clear_result = False, propagate_error = False):
        current_values = {"bkgrd": self.bkgrd, "Qmin": self.Qmin, "Qmax": self.Qmax}

        try:
            # Process Background, Qmin, and Qmax Inputs
            self.set_parameter("bkgrd", 'Background',self.bkgrd_input_box, 
                               default_val = 0, min_val = 0, 
                               error_msg = "Invalid background value. Must be >= 0.")
            self.set_parameter("Qmin", 'Q min', self.Qmin_input_box,
                               default_val = 0, min_val=0, 
                               error_msg = "Invalid Q min value. Must be >= 0.")
            self.set_parameter("Qmax", 'Q max', self.Qmax_input_box, 
                               default_val = np.inf, min_val=0, 
                               error_msg = "Invalid Q max value. Must be > 0.")
        except ValueError as e:
            if not propagate_error:
                self.show_error_message(str(e))
                return
            else:
                raise ValueError(str(e))

        # Replot if any values changed
        if any(current_values[attr] != getattr(self, attr) 
               for attr in current_values) or clear_result:            
            print(f"Background value set to: {self.bkgrd}")
            print(f"Minimum Q value set to: {self.Qmin}")
            print(f"Maximum Q value set to: {self.Qmax}")

            self.QQ_subtract, self.IQ_subtract, self.QQ_trim, self.IQ_trim = \
                bf.subtract_background(self.QQ_origin, self.IQ_origin, self.bkgrd, self.Qmin, self.Qmax)

            bf.plot_SANS_subtract(self.QQ_trim, self.IQ_trim, self.QQ_origin, self.bkgrd, self.figure_SAS, self.canvas_SAS)
            bf.plot_SANS_fit(self.QQ_trim, self.IQ_trim, self.figure_SAS_fitted, self.canvas_SAS_fitted, which='input')
            self.clear_result()
            self.recalc_PDSP_input_button.setEnabled(False)
            self.run_fit_button.setEnabled(True)
            self.save_result_button.setEnabled(False)
            
    # Set PDSP fit parameters: contrast, bulk density, pore radius for SSA extrapolate
    # number of points for SSA extrapolation
    def set_PDSP_fit_inputs(self, propagate_error = False):
        try:
            self.set_parameter('contrast', 'Contrast between 2 phases', 
                               self.contrast_input_box, default_val = 3e10, min_val = 0, 
                               error_msg = 'Invalid contrast value. Must be > 0')
            self.set_parameter('density', 'Sample bulk density', 
                               self.density_input_box, default_val = 1.0, min_val = 0, 
                               error_msg = 'Invalid density value. Must be > 0')
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
                                            'Must be > 0'))
            if abs(prev_val - self.r_SSA_extrapolate) < 1e-8:
                self.result_grid.removeWidget(self.SSA_label)
                self.SSA_label = QtWdgt.QLabel('SSA extrapolated<br>to {:.2f} nm'
                                               .format(self.r_SSA_extrapolate) +
                                               ' (cm<sup>2</sup>/cm<sup>3</sup>)')
                self.result_grid.addWidget(self.SSA_label, 3, 0)
        except ValueError as e:
            if not propagate_error:
                self.show_error_message(str(e))
                return
            else:
                raise ValueError(str(e))

    # Function used for determining whether the system majorly consists of void
    # or solid for porosity calculation. The input is read based on the state
    # of the flip button LeftRightSwitch()
    def set_major_phase(self, state):
        # Update the label when the switch toggles
        if state == 'Left':
            self.major_phase = 'solid'
        else:
            self.major_phase = 'void'
        print(f"Major phase of sample: {self.major_phase}")
     
    # This function fits, returns, and plots the PDSP result based on the provided
    # inputs. Activated when the button 'Fit PDSP Model!' is clicked.
    def run_fit_func(self):
        try:
            # Confirming fitting inputs
            self.set_bkgrd_Q_range(clear_result=True, propagate_error = True)
            self.set_PDSP_fit_inputs(propagate_error = True)
        except ValueError as e:
            self.show_error_message(str(e))
            return
            
        # Update program to show old results have been cleared
        QtWdgt.QApplication.processEvents()
        time.sleep(0.001)  
        
        # Run fit function
        try:
            self.rr, self.IQ_fitted, self.IQ0_fitted, self.f_r, self.f_dash_r, self.SSA,\
                self.dV_dr, self.phi, self.Vpore_avg, self.phi_on_Vavg, self.SSA_extrapolate = \
                    bf.fit_PDSP_model(self.QQ_trim, self.IQ_trim, self.pts_per_dec, 
                                      self.contrast, self.density, self.r_SSA_extrapolate, 
                                      self.num_pts_SSA_extrapolate, self.major_phase)
        except ValueError as e:
            self.show_error_message(str(e))
            return
        print('Done PDSP fit!')
        
        # Display result, enable/disable corresponding buttons to prevent accidental inputs.
        self.display_result()        
        self.recalc_PDSP_input_button.setEnabled(True)
        self.run_fit_button.setEnabled(False)
        self.save_result_button.setEnabled(True)
    
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
        self.porosity_box.setText('{:.3f}'.format(self.phi))
        self.pore_concentration_box.setText('{:.3e}'.format(self.phi_on_Vavg))
        self.pore_volume_box.setText('{:.3e}'.format(self.Vpore_avg))
        self.SSA_box.setText('{:.3e}'.format(self.SSA_extrapolate[0]))
        
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
        data_table = np.column_stack((self.rr, self.f_r, self.SSA, self.dV_dr))
        
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
                file.write('Contrast between 2 phases (cm-2): ' + '{:.3e}'.format(self.contrast))
                file.write('\n')
                file.write('Density of Solid (g/cm3): ' + '{:.3f}'.format(self.density))
                file.write('\n\n')            
                file.write('Porosity: ' + '{:.5e}'.format(self.phi))
                file.write('\n')
                file.write('Average Pore Volume (cm3): ' + '{:.5e}'.format(self.Vpore_avg))
                file.write('\n')
                file.write('Pore Concentration (cm-3): ' + '{:.5e}'.format(self.phi_on_Vavg))
                file.write('\n')
                file.write('SSA interpolated to r = {:.2f} nm (cm2/cm3): '
                           .format(self.r_SSA_extrapolate) +
                           '{:.5e}'.format(self.SSA_extrapolate[0]))
                file.write('\n\n')
                
                # Write result table for r vs f(r), SSA(R), and dV/dr
                file.write('Pore size distribution table\n')
                file.write('\t\t'.join(['r\t', 'f(r)', 'SSA\t', 'dV/dr']) + '\n')
                [file.write('\t'.join(val if isinstance(val, str)
                                      else '{:.5e}'.format(val) 
                                      if val > 0 
                                      else '{:.4e}'.format(val) 
                                      for val in line) + '\n')
                 for line in data_table]
            print(f"File saved as {file_name_save} in {self.chosen_save_folder_dir}")

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
    def show_error_message(self, message):
        msg_box = QtWdgt.QMessageBox()
        msg_box.setIcon(QtWdgt.QMessageBox.Warning)
        msg_box.setWindowTitle("Input Error")
        msg_box.setText(message)
        msg_box.setStandardButtons(QtWdgt.QMessageBox.Ok)
        msg_box.exec_()        
   
       
# Flip switch used for selecting either solid or void as the major phase of the 
# porous system.
class LeftRightSwitch(QtWdgt.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.button_width = self.find_button_size(88*parent.scale)
        self.button_height = self.find_button_size(40*parent.scale)
        self.current_state = "Left"  # Initial state
        self.setFixedSize(self.button_width, self.button_height)  # Fixed size for the switch
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)  # Cursor changes to a hand pointer

    def paintEvent(self, event):
        # Override the paintEvent to draw the switch
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)

        # Draw the background
        if self.current_state == "Left":
            painter.setBrush(QtGui.QBrush(QtGui.QColor(40, 40, 40)))  # Dark grey "Left"
        else:
            painter.setBrush(QtGui.QBrush(QtGui.QColor(40, 100, 200)))  # Light blue for "Right"

        painter.setPen(QtCore.Qt.NoPen)
        painter.drawRoundedRect(0, 0, self.width(), self.height(), 
                                int(self.button_height/2), int(self.button_height/2))

        # Draw the handle
        if self.current_state == "Left":
            painter.setBrush(QtGui.QBrush(QtCore.Qt.white))
            painter.drawEllipse(int(self.button_height*0.1), int(self.button_height*0.1), 
                                int(self.button_height*0.8), int(self.button_height*0.8))
        else:
            painter.setBrush(QtGui.QBrush(QtCore.Qt.white))
            painter.drawEllipse(self.width() - int(self.button_height*0.9), int(self.button_height*0.1), 
                                int(self.button_height*0.8), int(self.button_height*0.8))

    def mousePressEvent(self, event):
        # Toggle the state when the switch is clicked
        if event.button() == QtCore.Qt.LeftButton:
            self.current_state = "Right" if self.current_state == "Left" else "Left"
            self.update()  # Repaint the switch
            self.window().set_major_phase(self.current_state)  # Notify parent widget
            
    def find_button_size(self, size_start):
        """
        Find the smallest integer greater than or equal to `size_start` such that
        multiplying it by both 0.8 and 0.1 results in integers.
        Handles non-integer inputs as well.
        """
        # Start with the ceiling of the input to ensure it's rounded up to the nearest integer
        size = round(size_start)
        
        while True:
            if (size * 0.8).is_integer() and (size * 0.1).is_integer():
                return size
            size += 1
            
            
if __name__ == "__main__":
    app = QtWdgt.QApplication(sys.argv)
    main_window = PRINSAS_App()
    main_window.show()
    sys.exit(app.exec_())