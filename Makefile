# Compilateur utilise
CC = /usr/lib64/openmpi/bin/mpicxx

# options en mode optimise
OPTIM_FLAG = -O3 -std=c++11 -Wall -Woverloaded-virtual -I Eigen/Eigen

# options en mode debug
DEBUG_FLAG = -g -std=c++11 -Wall -Woverloaded-virtual -I Eigen/Eigen

# executable produit
PROG = run

# fichier source a compiler
SRC = Main.cc Functions.cpp SolverCG.cpp Timestep.cpp Tools.cpp

# par defaut on compile en optimise
optim : $(SRC)
	$(CC) $(SRC) $(OPTIM_FLAG) -o $(PROG)

debug : $(SRC)
	$(CC) $(SRC) $(DEBUG_FLAG) -o $(PROG)
