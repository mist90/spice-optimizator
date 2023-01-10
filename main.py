#    MIT License
#  Copyright (c) 2022 Mikhail Tegin
#  michail3110@gmail.com
#  
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

import mos3
import grapher
import matplotlib.pyplot as plt
import configparser
import json
import sys

def GetOption(config, Section, Parameter, Default):
    if config.has_option(Section, Parameter):
        return config.get(Section, Parameter)
    else:
        return Default

# Getting file name from command line
if len(sys.argv) > 1:
    fileName = sys.argv[1]
else:
    fileName = 'default.json'
print('Open config file:' + fileName)

# Reading configuration file
config = configparser.ConfigParser()
config.read(fileName)
# Mandatory parameters
optimize = config.get("Parameters", "optimize")
Vgs_list = list(map(float, json.loads(config.get("OutputGraph", "Vgs"))))
Vds_max_output = float(config.get("OutputGraph", "Vds_max"))
Vds_transfer = float(config.get("TransferGraph", "Vds"))
Vgs_max_transfer = float(config.get("TransferGraph", "Vgs_max"))
Ymin = float(config.get("TransferGraph", "Ymin"))
Ymax = float(config.get("TransferGraph", "Ymax"))

# Setup model
TestModel = mos3.Transistor()
TestModel.VT0 = float(GetOption(config, "Parameters", "VT0", 1.0))
TestModel.NSUB = float(GetOption(config, "Parameters", "NSUB", 1.0E+15))
TestModel.NFS = float(GetOption(config, "Parameters", "NFS", 1.0E+9))
TestModel.Rd = float(GetOption(config, "Parameters", "Rd", 0.0))
TestModel.Rs = float(GetOption(config, "Parameters", "Rs", 0.0))
TestModel.KP = float(GetOption(config, "Parameters", "KP", TestModel.KP))
TestModel.U0 = float(GetOption(config, "Parameters", "U0", TestModel.U0))
TestModel.Weff = float(GetOption(config, "Parameters", "Weff", TestModel.Weff))
TestModel.Leff = float(GetOption(config, "Parameters", "Leff", TestModel.Leff))
TestModel.TOX = float(GetOption(config, "Parameters", "TOX", TestModel.TOX))
TestModel.VMAX = float(GetOption(config, "Parameters", "VMAX", TestModel.VMAX))
TestModel.THETA = float(GetOption(config, "Parameters", "THETA", TestModel.THETA))
TestModel.XJ = float(GetOption(config, "Parameters", "XJ", TestModel.XJ))
TestModel.PHI = float(GetOption(config, "Parameters", "PHI", TestModel.PHI))
TestModel.DELTA = float(GetOption(config, "Parameters", "DELTA", TestModel.DELTA))
TestModel.GAMMA = float(GetOption(config, "Parameters", "GAMMA", TestModel.GAMMA))
TestModel.ETA = float(GetOption(config, "Parameters", "ETA", TestModel.ETA))
TestModel.KAPPA = float(GetOption(config, "Parameters", "KAPPA", TestModel.KAPPA))
TestModel.update()

Vgs = list(map(float, json.loads(config.get("TransferPoints", "Vgs"))))
Vds = list(map(float, json.loads(config.get("TransferPoints", "Vds"))))
Id = list(map(float, json.loads(config.get("TransferPoints", "Id"))))
if len(Vgs) != len(Vds) or len(Vds) != len(Id):
    print("Not correct points is assigned in section TransferPoints: len(Vgs)=" + str(len(Vgs)) + " len(Vds)=" + str(len(Vds)) + " len(Id)=" + str(len(Id)))
    exit(-1)
pointsSubTh = list(zip(Vgs, Vds, Id))
Vgs = list(map(float, json.loads(config.get("OutputPoints", "Vgs"))))
Vds = list(map(float, json.loads(config.get("OutputPoints", "Vds"))))
Id = list(map(float, json.loads(config.get("OutputPoints", "Id"))))
if len(Vgs) != len(Vds) or len(Vds) != len(Id):
    print("Not correct points is assigned in section OutputPoints: len(Vgs)=" + str(len(Vgs)) + " len(Vds)=" + str(len(Vds)) + " len(Id)=" + str(len(Id)))
    exit(-1)
pointsTh = list(zip(Vgs, Vds, Id))

if optimize == "Yes":
    # Extract parameters
    optimizer = mos3.Optimizer(TestModel, pointsSubTh, pointsTh)
    optimizer.VT0min = float(GetOption(config, "Parameters", "VT0min", optimizer.VT0min))
    optimizer.VT0max = float(GetOption(config, "Parameters", "VT0max", optimizer.VT0max))
    optimizer.NFSmin = float(GetOption(config, "Parameters", "NFSmin", optimizer.NFSmin))
    optimizer.NFSmax = float(GetOption(config, "Parameters", "NFSmax", optimizer.NFSmax))
    optimizer.KPmin = float(GetOption(config, "Parameters", "KPmin", optimizer.KPmin))
    optimizer.KPmax = float(GetOption(config, "Parameters", "KPmax", optimizer.KPmax))
    optimizer.THETAmin = float(GetOption(config, "Parameters", "THETAmin", optimizer.THETAmin))
    optimizer.THETAmax = float(GetOption(config, "Parameters", "THETAmax", optimizer.THETAmax))
    optimizer.KAPPAmin = float(GetOption(config, "Parameters", "KAPPAmin", optimizer.KAPPAmin))
    optimizer.KAPPAmax = float(GetOption(config, "Parameters", "KAPPAmax", optimizer.KAPPAmax))
    if not optimizer.run():
        print("Not Successful")
    print("Residuals: ", optimizer.residuals)


# Printing of Spice NMOS-parameters
print("VT0=" + str(TestModel.VT0))
print("NFS=%E" % TestModel.NFS)
print("KP=" + str(TestModel.KP))
print("THETA=" + str(TestModel.THETA))
print("KAPPA=" + str(TestModel.KAPPA))

# Plot graphs
fig = plt.figure()
axOutput = fig.add_subplot(2, 2, 1)
axOutput.set_title("Output characteristics")
axOutput.set_xlabel("Vds, V")
axOutput.set_ylabel("Ids, A")
axOutput.grid(True)
if len(Id) != 0:
    axOutput.scatter([point[1] for point in pointsTh], [point[2] for point in pointsTh])
graphOutput = grapher.OutputGrapher(TestModel, axOutput)
graphOutput.Vgs_list = Vgs_list
graphOutput.Vds_max = Vds_max_output
graphOutput.update()
axOutput.legend()

axTransfer = fig.add_subplot(2, 2, 2)
axTransfer.set_xlabel("Vgs, V")
axTransfer.set_ylabel("Ids, A")
axTransfer.set_ylim([Ymin, Ymax])
axTransfer.set_yscale("log")
axTransfer.grid(True)
graphTransfer = grapher.TransferGrapher(TestModel, axTransfer)
graphTransfer.Vds = Vds_transfer
graphTransfer.Vgs_max = Vgs_max_transfer
axTransfer.set_title("Transfer characteristic, Vds = " + str(Vds_transfer) + " V")
if len(Id) != 0:
    axTransfer.scatter([point[0] for point in pointsSubTh], [point[2] for point in pointsSubTh])
graphTransfer.update()

axTransfer = fig.add_subplot(2, 2, 4)
axTransfer.set_xlabel("Vgs, V")
axTransfer.set_ylabel("Ids, A")
axTransfer.set_ylim([Ymin, Ymax])
axTransfer.grid(True)
graphTransfer = grapher.TransferGrapher(TestModel, axTransfer)
graphTransfer.Vds = Vds_transfer
graphTransfer.Vgs_max = Vgs_max_transfer
axTransfer.set_title("Transfer characteristic, Vds = " + str(Vds_transfer) + " V")
if len(Id) != 0:
    axTransfer.scatter([point[0] for point in pointsSubTh], [point[2] for point in pointsSubTh])
graphTransfer.update()

plt.show()
