#! /bin/bash
g++ -g -I `root-config --incdir` gamma_addback.cc  `root-config --libs`
g++ -g -o dress_tree.o -I `root-config --incdir` dress_tree.cc  `root-config --libs`