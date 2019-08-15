INCDIRS = -I/usr/openwin/include -I/usr/local/include
WARNINGS = -Wall
CFLAGS = -g $(WARNINGS) $(INCDIRS)

libDataStructs.a:  Array.h CLinkedList.h DoubleRecurse.cpp DoubleRecurse.h \
		KdTree.cpp KdTree.h \
		ShellSort.h Stack.h
	g++ -c $(CFLAGS) $(LIBDIRS) \
		DoubleRecurse.cpp KdTree.cpp
	ar -r libDataStructs.a DoubleRecurse.o KdTree.o

