# Enceladus Biomolecule Simulator
By M.F.L. Meersman.
In order to use this model, please first contact m.f.l.meersman@student.tudelft.nl to ask for my permission.
If you just want to have a look at the processed data, please check out file 'dataframe.csv'.
All the produced graphs can be viewed and downloaded as '.PNG' images. It is however recommended to download the '.html'
files since these allow you to interact with the figures such that you can find and look up specific molecules.
In order to process the data yourself, please uncomment '# main()' of line 336 in file main.py, then run this file.
Going to 'https://rcbs.org' leads you to the international protein databank. From here, you can copy the
respective 'pdb-code-names' and paste them into one of three 'xxxxx_code_list.csv' files to process those respective
molecules.
This model simulaties the time required for the disintegration of possible lifeforms on Enceladus.
Going to 'dataframe.csv'-file will lead you to an overview of all the simulated molecule structures.
You can open this file best in Excel or any spreadsheet viewer. The meaning of all columns is:
- Identifier: This indicates whether the molecule is DNA, fatty-acid or proten/amino-acid
- code: This is the universal naming of the molecule that is simulated. Typing in this code in Google can give you additional info about the molecule.
- Required Height [m]: This is the distance you need to be elevated above the surface in order to capture the molecule before it brakes down.
- id: binary indication for Identifier
- Area [Å2]: This is the geometric area of the molecule structure, i.e. the area that will capture radiation flux.
- Time required [s]: This is the time from exiting the vent untill the molecule breaks down.
- URL: the link that leads to all the information of the molecule.
- Molecular complexity [-]: expressed as a function of number of atoms and variations of atoms.
