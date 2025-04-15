LOTF
====


.. image:: lotf.jpg
   :width: 300px
   :align: center

|

**LOTF** is a C++/ROOT track reconstruction algorithm designed for the Straw Tube Tracker in PANDA. If you wish to use it, please get in contact with me. This code was meant to be included in the Pandaroot reconstruction software and therefore requires several dependencies which are not public. 

LOTF was published in the following paper:

- `Reconstructing charged-particle trajectories in the PANDA Straw Tube Tracker using the LOcal Track Finder (LOTF) algorithm <https://link.springer.com/article/10.1140/epja/s10050-023-01005-8>`_ 



Paper summary
-------------

LOTF is a fast and robust track reconstruction algorithm optimized for the PANDA Straw Tube Tracker. It uses:

- A graph-based representation of the STT geometry
- A three-phase reconstruction approach: **connect**, **fit**, and **merge**
- A system of **virtual nodes** to recover z-information from skewed layers
- A low computational complexity, achieving **68,000 hits/s** without parallelization

Result figures:
----------------


1. STT Geometry and Skewed Layers:

   .. raw:: html

      <div style="display: flex; justify-content: space-between;">
          <img src="Publication material/Figures/STTGridVirt_XYPlane.png" width="45%" />
          <img src="Publication material/Figures/STTGrid_ZYPlane.png" width="45%" />
      </div>

2. Simulation setup:

   .. raw:: html

      <div style="display: flex; justify-content: space-between;">
          <img src="Publication material/Figures/ntrackpsimAll.png" width="85%" />
      </div>


3. Track reconstruction performance comparison:

   .. raw:: html

      <div style="display: flex; justify-content: space-between;">
          <img src="Publication material/Figures/PandaQA3GevV2.png" width="40%" />
          <img src="Publication material/Figures/PandaQA15GevV2.png" width="40%" />
      </div>

4. Event-mixing drop off performance

   .. raw:: html

      <div style="display: flex; justify-content: space-between;">
          <img src="Publication material/Figures/PandaQAEvtMixV2.png" width="40%" />
      </div>

5. Time and recontruction rate

   .. raw:: html


      <div style="display: flex; justify-content: space-between;">
          <img src="Publication material/Figures/AllTime.png" width="45%" />
          <img src="Publication material/Figures/Allhitsrate.png" width="45%" />
      </div>

Notes:
------
ROOT files containing the analysis results are available on Zenodo.
- `Root files <https://zenodo.org/records/15223160>`_ 

AUTHOR
------

- Simon Gazagnes <sgsgazagnes@gmail.com>
