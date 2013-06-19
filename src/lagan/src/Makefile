CC = gcc $(CFLAGS)
CPP = g++ $(CFLAGS)
CFLAGS = -O3 # -Wall -W
TRGT_DIR = ..

all: ../anchors ../chaos ../order ../mlagan ../prolagan ../utils/bin2mf ../utils/bin2bl ../utils/cextract ../utils/cstat ../utils/contigorder ../utils/getbounds ../utils/getlength ../utils/getoverlap ../utils/rc ../utils/seqmerge ../utils/scorealign ../utils/scorecontigs ../utils/getcontigpos ../utils/fa2xfa ../utils/Glue ../utils/dotplot ../utils/overlay
	(cd glocal; $(MAKE))
clean: 
	rm -f *.o *~ utils/*~ mlagan.purify core
	(cd glocal; $(MAKE) clean)
../anchors: anchors.c skiplist.c
	$(CC) -o $(TRGT_DIR)/anchors anchors.c skiplist.c
../chaos: fchaos.c thrtrie.c skiplist.c global.c translate.c mempage.c filebuffer.c
	$(CC) -o $(TRGT_DIR)/chaos fchaos.c thrtrie.c skiplist.c global.c translate.c filebuffer.c -lm -DCHAOS__FLAG
../order: order.c diagmatrix.c filebuffer.c
	$(CC) -o $(TRGT_DIR)/order order.c diagmatrix.c filebuffer.c
../mlagan: mlagan.c diagmatrix.c multial.c skiplist.c filebuffer.c
	$(CC) -o $(TRGT_DIR)/mlagan mlagan.c multial.c diagmatrix.c skiplist.c filebuffer.c -lm -DMULTIAL__FLAG
../prolagan: prolagan.c diagmatrix.c multial.c skiplist.c filebuffer.c
	$(CC) -o $(TRGT_DIR)/prolagan prolagan.c multial.c diagmatrix.c skiplist.c filebuffer.c -lm -DMULTIAL__FLAG
../utils/bin2mf: utils/bin2mf.c
	$(CC) -o $(TRGT_DIR)/utils/bin2mf utils/bin2mf.c
../utils/bin2bl: utils/bin2bl.c
	$(CC) -o $(TRGT_DIR)/utils/bin2bl utils/bin2bl.c
../utils/cextract: utils/cextract.c
	$(CC) -o $(TRGT_DIR)/utils/cextract utils/cextract.c
../utils/cstat: utils/cstat.c
	$(CC) -o $(TRGT_DIR)/utils/cstat utils/cstat.c
../utils/contigorder: utils/contigorder.c
	$(CC) -o $(TRGT_DIR)/utils/contigorder utils/contigorder.c
../utils/getbounds: utils/getbounds.c
	$(CC) -o $(TRGT_DIR)/utils/getbounds utils/getbounds.c
../utils/getcontigpos: utils/getcontigpos.c
	$(CC) -o $(TRGT_DIR)/utils/getcontigpos utils/getcontigpos.c
../utils/getlength: utils/getlength.c
	$(CC) -o $(TRGT_DIR)/utils/getlength utils/getlength.c
../utils/getoverlap: utils/getoverlap.c
	$(CC) -o $(TRGT_DIR)/utils/getoverlap utils/getoverlap.c
../utils/rc: utils/rc.c
	$(CC) -o $(TRGT_DIR)/utils/rc utils/rc.c
../utils/seqmerge: utils/seqmerge.c
	$(CC) -o $(TRGT_DIR)/utils/seqmerge utils/seqmerge.c
../utils/scorealign: utils/scorealign.c
	$(CC) -o $(TRGT_DIR)/utils/scorealign utils/scorealign.c -lm
../utils/scorecontigs: utils/scorecontigs.c
	$(CC) -o $(TRGT_DIR)/utils/scorecontigs utils/scorecontigs.c -lm
../utils/fa2xfa: utils/fa2xfa.c
	$(CC) -o $(TRGT_DIR)/utils/fa2xfa utils/fa2xfa.c
../utils/overlay: utils/overlay.c
	$(CC) -o $(TRGT_DIR)/utils/overlay utils/overlay.c
../utils/Glue: utils/Glue.cpp
	$(CPP) -o $(TRGT_DIR)/utils/Glue utils/Glue.cpp
../utils/dotplot: utils/dotplot.cpp
	$(CPP) -o $(TRGT_DIR)/utils/dotplot utils/dotplot.cpp
