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

import matplotlib.pyplot as plt

class OutputGrapher():
    def __init__(self, model, axis):
        self.numPoints = 200
        self.Vds_max = 25.0
        self.Vgs_list = [2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        self.model = model
        self.axis = axis
    def update(self):
        Vds_list = [self.Vds_max/float(self.numPoints)*float(i) for i in range(self.numPoints + 1)]
        for Vgs in self.Vgs_list:
            Id_list = [self.model.Id(Vgs, Vds) for Vds in Vds_list]
            self.axis.plot(Vds_list, Id_list, label = 'Vgs=' + str(Vgs) + ' V')

class TransferGrapher():
    def __init__(self, model, axis):
        self.numPoints = 200
        self.Vds = 50.0
        self.Vgs_max = 10.0
        self.model = model
        self.axis = axis
    def update(self):
        Vgs_list = [self.Vgs_max/float(self.numPoints)*float(i) for i in range(self.numPoints + 1)]
        Id_list = [self.model.Id(Vgs, self.Vds) for Vgs in Vgs_list]
        self.axis.plot(Vgs_list, Id_list)

