
from cc3d import CompuCellSetup
        

from CD8TcellProjectSteppables import CD8TcellProjectSteppable

CompuCellSetup.register_steppable(steppable=CD8TcellProjectSteppable(frequency=1))

from CD8TcellProjectSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

CompuCellSetup.run()
