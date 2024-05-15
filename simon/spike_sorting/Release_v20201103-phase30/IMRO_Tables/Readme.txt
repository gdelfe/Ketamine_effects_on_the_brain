------------------------------------------------------------
Premade pattern files for Neuropixels 2.0 single shank probe
------------------------------------------------------------

Single bank files – selects all sites in the bank specified in the name:
NPtype21_bank0_ref0.imro
NPtype21_bank1_ref0.imro
NPtype21_bank2_ref0.imro
NPtype21_bank0_ref1.imro    (uses tip reference)

Block of 384 channels starting from row specified in name:
NPtype21_botRow256.imro
NPtype21_botRow448.imro

Checkerboard pattern covering two banks, specified in name:
NPtype21_bank01_checker_ref0.imro


----------------------------------------------------------
Premade pattern files for Neuropixels 2.0 four shank probe
----------------------------------------------------------

Single shank files – selects 384 sites on one shank, starting from botRow:
NPtype24_shank0_botRow0_ref0.imro
NPtype24_shank0_botRow0_ref1.imro	(uses shank 0 tip reference)
NPtype24_shank2_botRow0_ref0.imro
NPtype24_shank2_botRow0_ref3.imro	(uses shank 2 tip reference)

Horizontal stripe files – 96-site-tall stripe, starting from botRow:
NPtype24_hStripe_botRow0_ref0.imro
NPtype24_hStripe_botRow272_ref0.imro
NPtype24_hStripe_botRow592_ref0.imro


-----
Notes
-----

To create variants of these patterns, we’ve provided two MATLAB scripts:

nptype21_imro.m
nptype24_imro.m

To run, edit the script to specify the pattern type and relevant parameters (described in the comments). When running SpikeGLX, use the Shank Viewer to check that the pattern matches your expectations.
