#!/usr/bin/perl

use warnings;
use Getopt::Long;


my %Opt;
if(!GetOptions(
	\%Opt, 
	'--h',
)){
	&print_usage();
	exit(1);
}

#if '-h',the help message will be printed
if($Opt{'h'}){
	&print_usage();
	 exit(1);
}

our ($input_msisensor_file,$output_vcf_file)=@ARGV;

if ( -e $input_msisensor_file ) {
	my $dis_pos=index($input_msisensor_file,"_dis");
	#my $slash_pos=rindex($input_msisensor_file,"/");
	our $soamtic_dis_file="";
	#our $inter_name= ""; 
	#our $mediate_vcf_output="";
	if($dis_pos > 0){
	#	$inter_name=substr($input_msisensor_file,$slash_pos+1,$dis_pos-$slash_pos-1);
	#	$mediate_vcf_output="MSISensor1_somatic_".$inter_name.".vcf";
	#	$soamtic_dis_file="MSISensor1_somatic_".$inter_name."_dis";
	#}elsif($dis_pos){
	#	$inter_name=substr($input_msisensor_file,0,$dis_pos);
	#	$mediate_vcf_output="MSISensor1_somatic_".$inter_name.".vcf";
	#	$soamtic_dis_file="MSISensor1_somatic_".$inter_name."_dis";
	}else{
		print "the input file must be the distribution from MSISensor1 and ended with '_dis' or '_dis.gz' \n";
		exit(1);
	}

	my $file_handle =  $input_msisensor_file =~ /gz$/ ? qq{gzip -dc $input_msisensor_file | } : qq{$input_msisensor_file};
	
	## convert somatic msi distribution file to somatic vcf file
	open(INPUT, $input_msisensor_file) or die $!;
	my @lines = <INPUT>;
	chomp(@lines);
	my $input_len=$#lines;
	open DELETE_SIETS, ">deleted_msi_sites" or die$!;
	open OUTPUT_MED, ">$output_vcf_file" or die $!;
	print OUTPUT_MED "##fileformat=VCFv4.2\n";
	print OUTPUT_MED "##INFO=<ID=CONEXTE_LEFT,Number=1,Type=String,Description=\"Five bases on  the left of MSI site\">\n";
	print OUTPUT_MED "##INFO=<ID=CONEXTE_RIGHT,Number=1,Type=String,Description=\"Five bases on  the right of MSI site\">\n";
	print OUTPUT_MED "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print OUTPUT_MED "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n";
	print OUTPUT_MED "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the alt alleles in the order listed by ALT\">\n";
	print OUTPUT_MED "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n";
	
	my $n=0;
	while($n <= $input_len){
		if($n %3 == 0){ 
                        my @ref_line=split(/ /,$lines[$n]);
			chomp(@ref_line);
                        our $chr=$ref_line[0];
                        our $pos=$ref_line[1];
                        our $context_left=$ref_line[2];
			our $context_right=$ref_line[4];
			our $repeat= $ref_line[3];
                        my $p=index($repeat,"[");
                        our $repeat_num=substr($repeat,0,$p);
                        my $o=index($repeat,"]");
                        our $repeat_base=substr($repeat,$p+1,$o-$p-1);
			our $ref=$repeat_base x $repeat_num;
                }elsif($n %3 ==1){
                        our @normal_line=split(/ /,$lines[$n]);
			chomp(@normal_line);
                        our $normal_depth=0;
			my $p=1; 
			while($p < @normal_line){
				my $tmp=$normal_line[$p];
				$normal_depth+=$tmp;
				$p++;
			}
		}else{
			my @tumor_line=split(/ /,$lines[$n]);
			chomp(@tumor_line);
			my $tumor_depth = 0;
			my $p=1;
			our %depth_pos_map=(); 
			while($p < @tumor_line){
				my $tmp=$tumor_line[$p];
				$tumor_depth += $tmp;
				# the dic of depth and letgh
				$depth_pos_map{$tmp} = $p;
				#print "map is : $tmp -> $depth_pos_map{$tmp}\n";
				$p++;
			}
			shift(@tumor_line);
			# sort the tumor array in descending order
			my @tumor_line_sorted = sort { $b <=> $a } @tumor_line;
			my $tumor_alt_dep = 0; # the depth of second repeated unit times
			my $tumor_repeatunit_length = 0 ; # the second max repeat unit length;
			if($depth_pos_map{$tumor_line_sorted[0]} eq $repeat_num){ # The max repeat num is equal to the reference repeat num 
				if($tumor_line_sorted[1] ne 0){ 
					# and the second max repeat num is not equal to 0 
					$tumor_alt_dep = $tumor_line_sorted[1];
					$tumor_repeatunit_length = $depth_pos_map{$tumor_line_sorted[1]};
				}else{
					# and the second max repeat num is equal to 0
					# $tumor_alt_dep = $tumor_line_sorted[0];
					$m=$n-2;
					$o=$n-1;
					print DELETE_SIETS "Please note that MSI Site: $lines[$m] is deleted. Because the MS distribution in Tumor is the same to the reference. The distributions are:\n";
					print DELETE_SIETS "$lines[$o]\n$lines[$n]\n";
					$n++;
					next();
				}	
			}else{
				$tumor_alt_dep = $tumor_line_sorted[0];
				$tumor_repeatunit_length = $depth_pos_map{$tumor_line_sorted[0]};
			}

			our $tumor_alt = "";
			our $type = "";
			#print "hexy: second_repeatunit_length is : $tumor_repeatunit_length\n";
                        if($tumor_repeatunit_length < $repeat_num) {
				# deletion
				$type = "DEL";
				my $deletion_time= $repeat_num - $tumor_repeatunit_length;
				my $tmp_ref = $repeat_base x $deletion_time;
				my $left_base = substr($context_left,4,1);
				$ref = $left_base.$tmp_ref;
				$tumor_alt = $left_base;
			}elsif($tumor_repeatunit_length > $repeat_num){
				# insertion
				$type = "INS";
				my $insertion_time= $tumor_repeatunit_length - $repeat_num;
                                my $tmp_ref = $repeat_base x $insertion_time;
                                my $left_base = substr($context_left,4,1);
                                $ref = $left_base;
                                $tumor_alt = $left_base.$tmp_ref;
			}
			print OUTPUT_MED "$chr\t$pos\t.\t$ref\t$tumor_alt\t.\t.\tCONEXTE_LEFT=$context_left;CONTEXT_RIGHT=$context_right;TYPE=$type\tGT:DP:AD\t0\/1:$tumor_depth:$tumor_alt_dep\n";	
			$normal_depth=0; $tumor_depth=0;
		}
		$n++;
	} # end of while
	close INPUT;
	close DELETE_SIETS;
	close OUTPUT_MED;

}else{ # end of if
	exit(1);
}

sub print_usage{
	print "<< Usage: perl MSISensor2Vepvcf.pl msi_dis_file vep_annotated_output path_of_vep path_of_vep_database
	[options] =
	[-h]	prints this message
	[msi_dis_file] 	Required	 the distribution file from MSISensor1
	[output]	Required	vep annoated file in VCF 
	[vep]	Requied	path of vep
	[vep_db]	Required	path of database or cache of vep\n"
} 
