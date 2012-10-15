//***************************************************************************
//* Copyright (c) 2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file COPYING for details.
//****************************************************************************

#include "outputgenerator.h"

namespace SyntenyFinder
{
	const std::string OutputGenerator::circosTemplate_ = "\
karyotype = circos.sequences.txt\n\
chromosomes_units = 1000000\n\
\n\
<ideogram>\n\
	<spacing>\n\
		default = 0.01r\n\
	</spacing>\n\
	radius    = 0.9r\n\
	thickness = 50p\n\
	fill      = yes\n\
	show_label       = yes\n\
	label_font       = default \n\
	label_radius     = dims(image,radius) - 80p\n\
	label_size       = 80\n\
	label_parallel   = yes\n\
</ideogram>\n\
\n\
show_ticks       = yes\n\
show_tick_labels = yes\n\
\n\
<ticks>\n\
	radius           = 1r\n\
	color            = black\n\
	thickness        = 4p\n\
	multiplier       = 1e-3\n\
	format           = %d\n\
	<tick>\n\
		spacing        = 0.05u\n\
		size           = 20p\n\
		show_label     = yes\n\
		label_size     = 15p\n\
		label_offset   = 10p\n\
	</tick>\n\
</ticks>\n\
\n\
<colors>\n\
	chr1* = red\n\
	chr2* = green\n\
	chr3* = blue\n\
	chr4* = orange\n\
</colors>\n\
\n\
<links>\n\
	crest = 1\n\
	bezier_radius = 0.3r\n\
	bezier_radius_purity = 1.0\n\
	<link>\n\
		show = yes\n\
		ribbon = yes\n\
		stroke_color = vdgrey\n\
		stroke_thickness = 2\n\
		file 		  = circos.segdup.txt\n\
		radius        = 0.99r\n\
		color         = red_a4\n\
		thickness     = 15\n\
	</link>\n\
</links>\n\
\n\
<highlights>\n\
	fill_color = green\n\
	<highlight>\n\
		file       = circos.highlight.txt\n\
		ideogram	= yes\n\
		fill_color = blue_a3\n\
		stroke_color = black\n\
		stroke_thickness = 4\n\
	</highlight>\n\
</highlights>\n\
\n\
################################################################\n\
# The remaining content is standard and required. It is imported \n\
# from default files in the Circos distribution.\n\
#\n\
# These should be present in every Circos configuration file and\n\
# overridden as required. To see the content of these files, \n\
# look in etc/ in the Circos distribution.\n\
\n\
<image>\n\
# Included from Circos distribution.\n\
<<include etc/image.conf>>\n\
</image>\n\
\n\
# RGB/HSV color definitions, color lists, location of fonts, fill patterns.\n\
# Included from Circos distribution.\n\
<<include etc/colors_fonts_patterns.conf>>\n\
\n\
# Debugging, I/O an dother system parameters\n\
# Included from Circos distribution.\n\
<<include etc/housekeeping.conf>>\n";

}