# --------------------------------------------------------------------------#
#								pySPHinXsys									#
# --------------------------------------------------------------------------#
# pySPHinXsys provides a python binding for SPHinXsys library.				#
#																			#
# Copyright(c) 2017-2024 The authors and the authors' affiliations.         #
#                                                                           #
# Licensed under the Apache License, Version 2.0 (the "License"); you may   #
# not use this file except in compliance with the License. You may obtain a #
# copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        #
#                                                                           #
# --------------------------------------------------------------------------#
# @file 	dambreak.y
# @brief 	dam break test with python binding. 
# @author	Chi ZHANG
# @version  1.0 
#			Dam break flow. 
#			-- Chi ZHANG

import pysphinxsys as pySPH

def main():
    dambreak = pySPH.Dambreak()
    dambreak.runCase(20)

if __name__ == "__main__":
    main()