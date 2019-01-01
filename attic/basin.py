class Basin(object):

    def __init__(self, datafile):


        return self

    def build(self):

        t0 = time.clock()
        lib.logger.TIMER_LEVEL = 3



        if self.verbose >= logger.WARN:
            self.check_sanity()
        if self.verbose > logger.NOTE:
            self.dump_input()



        return self

    kernel = build

if __name__ == '__main__':
    name = 'h2o.chk'
    bas = Basin(name)
    bas.verbose = 4
    bas.nrad = 101
    bas.iqudr = 'legendre'
    bas.mapr = 'becke'
    bas.bnrad = 101
    bas.bnpang = 5810
    bas.biqudr = 'legendre'
    bas.bmapr = 'becke'
    bas.non0tab = False

    bas.inuc = 0
    bas.kernel()

    bas.inuc = 1
    bas.kernel()
 
    bas.inuc = 2
    bas.kernel()

