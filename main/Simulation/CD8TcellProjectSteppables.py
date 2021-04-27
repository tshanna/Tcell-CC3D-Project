
# TODO: Modularize plots
# TOFIX: Simulation will crash if cellOI dies during the simulation. 
#        cellOI is used for plotting purposes

from cc3d.core.PySteppables import *
import numpy as np
import random
import sys
import math

ODE_vars = ['IR','IRa','Tb','Fs','Fsa','C']

small_num = sys.float_info.min
sec_per_mcs = 60
pi = math.pi

class CD8TcellProjectSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
        self.pop_selected = False
        
        self.did_seeding = False
        self.did_seeding_again = False
        
        self.cellOI = None

    def start(self):
        
        self.plot_win = self.add_new_plot_window(title='Intracellular ODE Values',
                                                 x_axis_title='Hours',
                                                 y_axis_title='Amount', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False,config_options={'legend': True})
        

        # self.plot_win2 = self.add_new_plot_window(title='Activated Cells',
                                                 # x_axis_title='Days',
                                                 # y_axis_title='Amount', x_scale_type='linear', y_scale_type='linear',
                                                 # grid=False,config_options={'legend': True})
                                                 
                                                 
        # self.plot_win3 = self.add_new_plot_window(title='Effector Cells',
                                                 # x_axis_title='Days',
                                                 # y_axis_title='Amount', x_scale_type='linear', y_scale_type='linear',
                                                 # grid=False,config_options={'legend': True})                                              
        
        self.plot_win4 = self.add_new_plot_window(title='Cell Count',
                                                 x_axis_title='Days post infection',
                                                 y_axis_title='Number of Cells', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False,config_options={'legend': True}) 
     
        self.plot_win.add_plot("IRa", style='Dots', color='green', size=4)
        self.plot_win.add_plot("Casp", style='Dots', color='cyan', size=4)
        # self.plot_win.add_plot("IL2cm", style='Dots', color='purple', size=4)
        # self.plot_win.add_plot("fAPC", style='Dots', color='purple', size=4)
        # self.plot_win.add_plot("IR", style='Dots', color='purple', size=4)
        # self.plot_win.add_plot("Fsa", style='Dots', color='blue', size=4)
        # self.plot_win.add_plot("Fs", style='Dots', color='red', size=4)
        self.plot_win.add_plot("Tb", style='Dots', color='yellow', size=4)
        
        # self.plot_win2.add_plot("Tb", style='Dots', color='yellow', size=4)
        # self.plot_win2.add_plot("Casp", style='Dots', color='cyan', size=4)
        
        # self.plot_win3.add_plot("Casp", style='Dots', color='cyan', size=4)   
        
        self.plot_win4.add_plot(plot_name='N',style='Lines', color='cyan',size=4)
        self.plot_win4.add_plot(plot_name='APC',style='Lines', color='red',size=4)
        self.plot_win4.add_plot(plot_name='A',style='Lines', color='yellow',size=4)
        self.plot_win4.add_plot(plot_name='E',style='Lines', color='purple',size=4)
        self.plot_win4.add_plot(plot_name='P',style='Lines', color='green',size=4)
        self.plot_win4.add_plot(plot_name='T',style='Lines', color='blue',size=4)
        

    def add_steering_panel(self):
        
        self.add_steering_param(name='Start', val=0, min_val=0, max_val=1, widget_name='slider', decimal_precision=0)

        self.add_steering_param(name='Initial Cell Population', val=33, min_val=4, max_val=100, widget_name='slider')
        self.add_steering_param(name='Initial APC Population', val=3, min_val=1, max_val=10, widget_name='slider')
        self.add_steering_param(name='lamc1: caspase feedback strength', val=0.01, min_val=0.0, max_val=0.03, decimal_precision=3, widget_name='slider')
        self.add_steering_param(name='lamT3: Tbet feedback strength', val= 0.01, min_val=0.01, max_val=1.5, decimal_precision=2, widget_name='slider')
        
    def process_steering_panel_data(self):
        
        if not self.did_seeding:
            start = self.get_steering_param('Start')
            if start == 1:
                self.did_seeding = True
                pop_cells = self.get_steering_param('Initial Cell Population')
                apc_cells = self.get_steering_param('Initial APC Population')
                self.populate_cells(pop_cells, apc_cells)

        
        
        
        for cell in self.cell_list_by_type(self.NAIVE, self.PREACTIVATED, self.ACTIVATED, self.EFFECTOR):
            cell.sbml.dp['lamc1'] = self.get_steering_param('lamc1: caspase feedback strength')
            cell.sbml.dp['lamT3'] = self.get_steering_param('lamT3: Tbet feedback strength')

    def populate_cells(self, pop_num, apc_num):
        
        model_string = """
        
        E1: -> IR ; lamR1 * ( fAPC ) + ( lamR2 + muDIL2 ) * IRa - ( muAIL2 )*( IL2cm )* IR - ( kR * IR ) ;           // Non-activated IL2R
        E2: -> IRa ; ( muAIL2 )*( IL2cm )* IR - (muDIL2 * IRa) - ( kRa * IRa ) ;                                     // Activated IL2R
        E3: -> Tb ; lamT1 * ( fAPC ) + lamT2 * ( (Tb) / (lamT3 + Tb) ) * Tb - ( kT * Tb ) ;                          // Tbet
        E4: -> Fs ; lamF + ( muDF * Fsa ) - H * ( muAF * Tbcm * Fs ) - ( kF * Fs ) ;                                 // Nn-activated Fas
        E5: -> Fsa ; H * ( muAF * Tbcm * Fs ) - ( muDF * Fsa ) - ( kFa * Fsa ) ;                                     // Activated Fas
        E6: -> C ; lamc1 * ( 1/ (1 + (lamc2*IRa) ) ) * ( 1/ (1 + (lamc3*fAPC) ) ) + ( lamc4 * Fsa ) - ( kC * C ) ;  // Caspase
        
        // Decay rates
        
        kR = 0.0029 ; // 1/min
        kRa = 0.0029 ; // 1/min
        kT = 0.0035 ; // 1/min
        kF = 0.0047 ; // 1/min
        kFa = 0.0047 ; // 1/min
        kC = 0.0038 ; // 1/min
        
        // Feedback strengths
        
        lamR1 = 0.0158 ; // M 1/min
        lamR2 = 0.001 ; // 1/min
        lamT1 = 0.01 ; // M 1/min
        lamT2 = 0.004 ; // 1/min
        lamT3 = 0.01 ; // M
        lamc1 = 0.01 ; // M 1/min
        lamc2 = 100 ; // 1/M
        lamc3 = 0.01 ; // N/A
        lamc4 = 0.004 ; // 1/min
        
        // Association/Dissociation rates
        
        muAIL2 = 6E8 ; // 1/(M*min)
        muDIL2 = 0.006 ; // 1/min
        muAF = 0.86E5 ; // 1/(M*min)
        muDF = 0.006 ; // 1/min
        
        lamF = 3.47E-5 ; // M 1/min
        
        // Initial conditions
        
        IL2cm = 0
        fAPC = 0
        Tbcm = 0
        H = 0
        
        IR0 = 0
        IR = IR0
        IRa0 = 0
        IRa = IRa0
        Tb0 = 0
        Tb = Tb0
        Fs0 = 0
        Fs = Fs0
        Fsa0 = 0
        Fsa = Fsa0
        C0 = 0
        C = C0
        
        """
        lamR3 = 1E-12
        
        #initialize cells throughout the domain
        i = 0
        while i < pop_num:
            x1 = random.sample(range(0,201),1)
            x = x1[0]
            y1 = random.sample(range(0,201),1)
            y = y1[0]
            z = 1
            
            # make sure cells are not overwriting eachother
            if self.cell_field[x,y,z] is None:
                cell = self.new_cell(self.NAIVE)
                self.cell_field[x, y, z] = cell
                i += 1
        
        # L = random.sample(range(0, pop_num), 3)
        L = random.sample(range(0, pop_num), apc_num)
      
        
        for cell in self.cell_list:
            print(cell.id)

            # initialize the 3 starting APC 
            if cell.id in L:
                cell.type = self.APC
                
                if cell.type == self.APC:

                    cell.targetVolume = 250
                    cell.lambdaVolume = 10
                    
                    life = random.sample(range(32,41),1)
                    
                    time = life[0] * sec_per_mcs
                    
                    cell.dict['lifespan'] = time
                
            else:
                
                cell.targetVolume = 25
                cell.lambdaVolume = 10
                cell.dict['lamR3'] = lamR3
                
                self.add_antimony_to_cell(model_string=model_string,
                                  model_name='dp',
                                  cell=cell,
                                  step_size=1)
    
    def step(self,mcs):

        if mcs == 0:
            while not self.pop_selected:
                if self.steering_param_dirty():
                    break 
                pass

        if not self.cellOI:
            for cell in self.cell_list_by_type(self.PREACTIVATED):
                if cell.sbml.dp['fAPC'] > 0:
                    self.cellOI = cell
                    break   
            
        if not mcs % 10:
            for cell in self.cell_list:
                if self.cellOI:
                    self.plot_win.add_data_point("IRa", mcs/60, self.cellOI.sbml.dp['IRa'])
                    self.plot_win.add_data_point("Tb", mcs/60, self.cellOI.sbml.dp['Tb']) 
                    self.plot_win.add_data_point("Casp", mcs/60, self.cellOI.sbml.dp['C'])
                    
                    # self.plot_win.add_data_point("fAPC", mcs, self.cellOI.sbml.dp['fAPC'])
            # for cell in self.cell_list_by_type(self.ACTIVATED):
                # if self.cellOI:
                    # self.plot_win2.add_data_point("Tb", mcs, self.cellOI.sbml.dp['Tb']) 
                    # self.plot_win2.add_data_point("Casp", mcs, self.cellOI.sbml.dp['C'])
            # for cell in self.cell_list_by_type(self.EFFECTOR):
                # if self.cellOI:
                    # self.plot_win3.add_data_point("Casp", mcs, self.cellOI.sbml.dp['C'])
        
        # if not self.cellOI:
            # for cell in self.cell_list_by_type(self.PREACTIVATED):
                # if cell.sbml.dp['fAPC'] > 0:
                    # self.cellOI = cell.id
                    # break
        
        # if self.cellOI:
            # if not mcs % 50 :
                # cellOI = self.fetch_cell_by_id(self.cellOI)
                # if cellOI is not None:
                    # self.plot_win.add_data_point("IRa", mcs, cellOI.sbml.dp['IRa'])
                     #self.plot_win.add_data_point("IR", mcs, cellOI.sbml.dp['IR'])
                    
                    # self.plot_win.add_data_point("Casp", mcs, cellOI.sbml.dp['C'])
                    # self.plot_win.add_data_point("Fs", mcs, cellOI.sbml.dp['Fs'])
                    # self.plot_win.add_data_point("Fsa", mcs, cellOI.sbml.dp['Fsa'])
                    # #self.plot_win.add_data_point("Il2cm", mcs, cellOI.sbml.dp['IL2cm'])
                    # self.plot_win.add_data_point("fAPC", mcs, cellOI.sbml.dp['fAPC'])
            
                    # self.plot_win2.add_data_point("Tb", mcs, cellOI.sbml.dp['Tb']) 
                # else:
                    # self.cellOI = None
        
        if not mcs % 10 :
            self.plot_win4.add_data_point("N", (mcs/1440)+3, len(self.cell_list_by_type(self.NAIVE)))
            self.plot_win4.add_data_point("P", (mcs/1440)+3, len(self.cell_list_by_type(self.PREACTIVATED)))
            self.plot_win4.add_data_point("A", (mcs/1440)+3, len(self.cell_list_by_type(self.ACTIVATED)))
            self.plot_win4.add_data_point("E", (mcs/1440)+3, len(self.cell_list_by_type(self.EFFECTOR)))
            self.plot_win4.add_data_point("APC", (mcs/1440)+3, len(self.cell_list_by_type(self.APC)))
            self.plot_win4.add_data_point("T", (mcs/1440)+3, len(self.cell_list))
            
        IL2_secretor = self.get_field_secretor('IL2')
        
        IL2 = self.field.IL2
        
        #lamR3 = 1E-12
        lamR4 = 0.0
        lam1 = 1E-12
        lamT4 = 0.0
        
        for cell in self.cell_list_by_type(self.APC):
            
            # movement for APC updated every 90 minutes (mcs)
            if mcs % 90 == 0:
                n = random.uniform(-pi,pi)
                # randomize if movement vector is positive or negative
                sign = random.uniform(-3,3)
                abs_sign = abs(sign)
                
                # ensure there is no divide by 0
                if sign == 0:
                    sign = sign + 1
                    abs_sign = abs(sign)
                
                x_temp = math.cos(n)
                y_temp = math.sin(n)
                
                # r = 20 for APC
                x = 20.0*(x_temp*x_temp)*(sign/abs_sign)
                y = 20.0*(y_temp*y_temp)*(sign/abs_sign)

                cell.lambdaVecX = x
                cell.lambdaVecY = y
            
            # APC death after lifespan is reached
            if mcs >= cell.dict['lifespan']:
                
                cell.targetVolume = 0.0
        
        for cell in self.cell_list_by_type(self.NAIVE, self.PREACTIVATED, self.ACTIVATED, self.EFFECTOR):
            
            if cell.targetVolume > 0.0:
                
                # movement for T cells updated every 90 minutes (mcs)
                if mcs % 90 == 0:
                    n = random.uniform(-pi,pi)
                    
                    # randomize if movement vector is positive or negative
                    sign = random.uniform(-3,3)
                    abs_sign = abs(sign)
                    
                    # ensure there is no divide by 0
                    if sign == 0:
                        sign = sign + 1
                        abs_sign = abs(sign)
                    
                    x_temp = math.cos(n)
                    y_temp = math.sin(n)
                    
                    x = 150.0*(x_temp*x_temp)*(sign/abs_sign)
                    y = 150.0*(y_temp*y_temp)*(sign/abs_sign)
                    
                    cell.lambdaVecX = x
                    cell.lambdaVecY = y            

                fAPC = 0
                Tbcm = 0
                cell.sbml.dp['H'] = 0
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        # count amount of APC a T cell is in contact with
                        if neighbor.type == self.APC:
                            fAPC += 1
                        Tbcm += cell.sbml.dp['Tbcm']
                        # track effector-effector or effector-activated contacts (heaviside function H)
                        if cell.type == self.EFFECTOR:
                            if neighbor.type == self.EFFECTOR or self.ACTIVATED:
                                cell.sbml.dp['H'] = 1
          
                il2_cm = IL2_secretor.amountSeenByCell(cell)

                #update ODE
                cell.sbml.dp['IL2cm'] = il2_cm
                cell.sbml.dp['fAPC'] = fAPC
                cell.sbml.dp['Tbcm'] = Tbcm    
                
                # print("lamc1 VALUE IS NOW: ", cell.sbml.dp['lamc1'])
                # print("lamT3 VALUE IS NOW: ", cell.sbml.dp['lamT3'])
                
                # second term PDE
                #secrete = ( cell.dict['lamR3']*(( cell.sbml.dp['IRa'] )/( lamR4 + cell.sbml.dp['IRa'] + small_num)) + lam1 * cell.sbml.dp['fAPC'] ) * ( 1 / (1 + lamT4 * cell.sbml.dp['Tb']) )/(cell.surface)
                secrete = lam1*cell.sbml.dp['fAPC']
                
                if cell.type is not self.NAIVE:
                    # secretion of IL2 by T cells 
                    IL2_secretor.secreteOutsideCellAtBoundary(cell, secrete)
                    
                    #Caspase threshold for effector, activated, preactivated
                    if cell.sbml.dp['C'] > 2.63:
                       cell.targetVolume = 0.0

                # if Tbet threshold reached, A -> E
                if cell.type == self.ACTIVATED:
                    if cell.sbml.dp['Tb'] > 40:
                        cell.type = self.EFFECTOR
                
                # Preactivated cells stop moving until activated
                if cell.type == self.PREACTIVATED:

                    cell.lambdaVecX = 0.0
                    cell.lambdaVecY = 0.0
                    
                    #if IL2 threshold reached, PA -> A
                    if cell.sbml.dp['IRa'] > 7: 
                        cell.type = self.ACTIVATED                 
                        #cell.dict['lamR3'] = 0.0
                
                # Naive -> PA
                if cell.type == self.NAIVE and fAPC > 0:
                    cell.type = self.PREACTIVATED
        
        # step the simulation
        self.timestep_sbml()            
       
    def finish(self):
        
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return



class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list_by_type(self.EFFECTOR, self.ACTIVATED):
            # Effector and activated t cells divide every ~8 hours
            if mcs % 480 == 0:
                cells_to_divide.append(cell)
        
        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        
        # one cell inherits k = [0.7,1.0] fraction of parent cell
        # other cell has 2 - k 
        
        # self.parent_cell.targetVolume = self.cell.volume * k                 

        self.clone_parent_2_child()    

        for v in ODE_vars:
            k = random.uniform(0.7,1.0)
            x = self.parent_cell.sbml.dp[v]
            self.child_cell.sbml.dp[v] = x*(2-k)
            self.parent_cell.sbml.dp[v] = x*k

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        #if self.parent_cell.type==1:
        #    self.child_cell.type=2
        #else:
        #    self.child_cell.type=1
        



