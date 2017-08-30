from GenericElement import GenericElement
from WaveguideElement import WaveguideElement
from WaveguideJunction import WaveguideJunction
import Utils
from matplotlib import pyplot as plt

class Filter(GenericElement):
    """
    Describes waveguide filter object

    Waveguide filter consists of ``WaveguideElement`` instances connected via
    ``WaveguideJunction`` elements. 
    """

    s11 = None
    s12 = None
    s21 = None
    s22 = None
    
    def __init__(self, descr, minModes, frequency):
        """
        Initialize Filter object

        descr     : a list of tuples (width, length) (geometric parameters of ``WaveguideElement``)
        minModes  : number of modes used in a most narrow section of the filter. 
        frequency : operating frequency
        """

        minWidth = min(map(lambda x: Utils.toSI(x[0]), descr))
        loc_descr = []
        for d in descr:
            width, length = Utils.toSI(d[0], d[1])
            numModes = int(minModes * width/minWidth)
            loc_descr.append((width, length, numModes))

        self.descr = []
        for i in xrange(len(loc_descr)):
            widthA, lengthA, numModesA = loc_descr[i]
            self.descr.append(("waveguide", widthA, lengthA, numModesA))
            try:
                widthB, lengthB, numModesB = loc_descr[i+1]
                self.descr.append(("junction", widthA, widthB, numModesA, numModesB))
            except IndexError:
                pass

        self.frequency = Utils.toSI(frequency)
        self.initialize_networks()
        self.update(self.frequency)

    def __repr__(self):
        fmt0 = "Filter object, consists of following waveguide sections:"
        fmt1 = "    {i}. width = {width}m, length = {length}m, num. modes = {numModes}"
        message = []
        message.append(fmt0)
        i = 0
        for item in self.descr:
            if item[0] == 'waveguide':
                i += 1
                width, length, numModes = item[1:]
                message.append(fmt1.format(i = i,
                                           width = width,
                                           length = length,
                                           numModes = numModes))
        return '\n'.join(message)

    def initialize_networks(self):
        self.networks = []
        for descr in self.descr:
            wtype = descr[0]
            if wtype == 'waveguide':
                width, length, numModes = descr[1:]
                network = WaveguideElement(width, length, numModes, self.frequency)
            elif wtype == 'junction':
                widthA, widthB, numModesA, numModesB = descr[1:]
                network = WaveguideJunction(widthA, widthB, numModesA, numModesB, self.frequency)
            else:
                pass

            self.networks.append(network)

    def update(self, frequency):
        self.frequency = Utils.toSI(frequency)
        result = None
        for network in self.networks:
            network.update(frequency)
            if not result:
                result = network
            else:
                result = result * network

        self.s11 = result.s11
        self.s12 = result.s12
        self.s21 = result.s21
        self.s22 = result.s22
                

def main():
    """
    Run example code (a.k.a. Filter6_30a)
    """
    
    import numpy as np
    #from matplotlib import pyplot as plt
    #from progressbar import ProgressBar
    import sys
    from time import time
    
    a	= "8.636mm"
    b	= "3.556mm"
    c1	= "3.375mm"
    c2	= "3.563mm"
    c3	= "3.784mm"
    c4	= "4.924mm"
    d1	= "4.925mm"
    d2	= "4.722mm"
    d3	= "3.77mm"
    dd	= "2.5mm"
    a1	= "7.112mm"
    ddd	= "6mm"
    descr = [
        (a1, ddd),
        (c4, dd),
        (a, d3),
        (c3, dd),
        (a, d2),
        (c2, dd),
        (a, d1),
        (c1, dd),
        (a, d1),
        (c2, dd),
        (a, d2),
        (c3, dd),
        (a, d3),
        (c4, dd),
        (a1, ddd)
        ]
    descr_SI = []
    for width, length in descr:
        descr_SI.append(Utils.toSI(width, length))

    #plt.figure()
    freqs = np.linspace(26e9, 36e9, num=101)
    nmodes = 10

    benchmark = False
    
    def calc_filter(descr_SI, nmodes, freqs):
        s21 = np.zeros(freqs.shape)
        #bar = ProgressBar(maxval = len(freqs))
        #bar.start()
        first = True
        f = Filter(descr_SI, nmodes, freqs[0])
        for i, freq in enumerate(freqs):
            f.update(freq)
            s21[i] = 20*np.log10(np.abs(f.s21[0,0]))
            #bar.update(i)
        #bar.finish()

        if not benchmark:
            plt.plot(freqs * 1e-9, s21, label = ("Num.modes = %i" % nmodes))

    calc_filter(descr_SI, nmodes, freqs)

    if not benchmark:
        plt.legend()
        plt.grid()
        plt.show()


if __name__ == "__main__":
    main()
