Given a set of D documents, with each document classified according to S subjects, and a number C
 of cabinets, assign documents to cabinets according to the similarity of subjects.
 The classification of a document is obtained by giving it a score, a positive number, to each subject,
 thus obtaining for each document a vector with S dimensions.
 The decision of which document goes into which cabinet is made by trying to minimize the overall
 “distances” of the subjects of the document assigned to a cabinet (maximize subject similarity).

Program takes up to 3 arguments: <input file> <output file> [cabinet]

OpenMp, openMPI , C are used in this project
