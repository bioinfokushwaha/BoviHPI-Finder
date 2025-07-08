#!/usr/bin/perl -w
use strict; use warnings;
use English; use Getopt::Long; use Pod::Usage; use File::Basename;

=head1 NAME

RGPred

=head1 SYNOPSIS


		                  ###################################################################################################################
                    ######################## *** RGPred: Resistance gene prediction and sequence extraction pipeline  *** ######################
                    		  ###################################################################################################################

		  perl RGPred.pl --Infile=FVXD.fna --Seqtype=T --Score_cutof=0.1 --Tempfile=YES 
        	OR
		  ./RGPred.pl --Infile=Protein.Pep --Seqtype=P --Score_cutof=0.1 --Tempfile=NO
		OR
		  ./RGPred.pl --Infile=genomic_seq.fasta --Seqtype=G --Score_cutof=0.1 --Tempfile=NO --species=arabidopsis

		Input Options:
		--Infile	:Genome/Transcriptome/Protein fasta sequences
		--Seqtype	:Type of sequences i.e. G:(Genome),T:(Transcriptome), P:(Protein)	  
		--Cutoff	:Threshold for RGene prediction (0, 0.1, 0.2,....1....)
		--Species	:Name of specie for ab-initio gene prediction  
		--Tempfile 	:YES if you want keep intermediate files, otherwise NO
		--help		:For help
		--stdout	:Print to stdout instead of files
		--(no)verbose	:Be verbose; default no.

		Output Files:
		.CDS			:Extracted nucleotide sequences from genome and Transcriptome.
		.PEP			:Protein sequences of corresponding CDS sequences from genome and Transcriptome.
    		.RPredResult		:RGene prediction results as table.
    		.CDS_Ext_RGeneSeq	:Nucleotide sequences of predicted Rgene.
    		.PEP_Ext_RGeneSeq	:Protein sequences of predicted Rgene of corresponding nucleotide sequences.
    	 


		RGPred Requirments: Unix environment, perl, bioperl
		        Suggestion: 1. To avoid running troubles, run RGPred script from inside RGCap folder and keep sequence file at same place 
			            2. Split your sequences into multiple sequence files if it is more than 40K


=head1 AUTHOR

sandeep.kushwaha@slu.se

=cut

my $help = 0; my $Infile = 0; my $Seqtype=0; my @head =''; my $stdout = 0; my $verbose = 0; my @probe_array=0; 
my $Tempfile = ''; my $Pep_file=''; my $cds_file=''; my $Basename=''; my $Species=''; my $Score_cutoff =''; 

GetOptions(
  		"Tempfile:s" => \$Tempfile,
		"Infile:s" => \$Infile,
  		"help!" => \$help,
  		"Seqtype:s" => \$Seqtype,
		"Species:s" => \$Species,
		#"Protein" => \$Seqtype,
		"Score_cutoff=f" => \$Score_cutoff,
  		"stdout!" => \$stdout,
  		"verbose!" => \$verbose,
	);

pod2usage(0) if $help;
pod2usage(-msg => "Error: please check to make sure the file exists and is a readable plain text file.\n", -exitval => 1) unless $Infile;
pod2usage(-msg => "Need either a sequence type ", -exitval => 1) unless $Seqtype;
pod2usage(-msg => "Need either a cutoff score for R-Gene Prediction ", -exitval => 1) unless $Score_cutoff;
pod2usage(-msg => "Need a intermediate file option or specify --stdout", -exitval => 1) unless $Tempfile;

BEGIN { our $start_run = time(); }

	 my ($basename, $path, $suffix) = fileparse($Infile, '\.[^\.]*');
	 $Basename=$basename;  chomp $Basename; print "Input basefile name:$Basename\n"; 


#################********: RPred Pipeline processing for Proteom,Transcriptome and Genome sequence :********#########################
	       ###############################################################################################

	if ($Seqtype eq "T")  { #print "Input file name $Infile \n"; 
				&filter_by_len($Infile,$Basename);
				&filter_duplicate($Infile,$Basename);
				&Transdecoder($Basename,$Tempfile);

				my $Pep_file=$Basename.".transdecoder.pep"; #print "Transdecoder protein file:$Pep_file basefile:$Basename.PEP\n"; sleep (2);
				&Format_Ambiguity_removal($Pep_file,"T"); 
 				system("mv $Basename.transdecoder.pep.Ncbi_ambiguity_clean $Basename.PEP");

				my $Cds_file=$Basename.".transdecoder.cds"; #print "Transdecoder CDS file :$Cds_file basefile:$Basename.CDS\n";sleep (2);
				&Format_Ambiguity_removal($Cds_file, "T");
 				system("mv $Basename.transdecoder.cds.Ncbi_ambiguity_clean $Basename.CDS");

				&Rpred($Basename,$Score_cutoff,$Tempfile);
				&SeqExtaction_Trans($Basename,$Tempfile);
			      }

	elsif ($Seqtype eq "P") { 
				 #print "\n Input file name $Infile\n"; sleep (2);
				 &filter_by_len($Infile,$Basename);
				 &filter_duplicate($Infile,$Basename);

				 &Format_Ambiguity_removal($Basename,"P"); 
				 $Basename=$basename;  chomp $Basename;
				 system("mv $Basename.Ncbi_ambiguity_clean $Basename.PEP");
				
				 &Rpred($Basename); 
				 &SeqExtaction_Prot($Pep_file,$Tempfile); 
				}

	elsif ($Seqtype eq "G") { if ($Species eq "") {print "\nError:No input for species OR Not able to read input for species\nPlease check carefully and re-run program.\n";}
				  
				  &Augustus($Infile,$Basename, $Species); 
			  	  my $aug_pep=$Basename.".aa";	
			  	  &Format_Ambiguity_removal($aug_pep,$Basename,'G'); 
				  system("mv $Basename.aa.Ncbi_ambiguity_clean $Basename.PEP");
				
				  my $aug_cds=$Basename.".codingseq";	
				  &Format_Ambiguity_removal($aug_cds,$Basename,'G'); 
				  system("mv $Basename.codingseq.Ncbi_ambiguity_clean $Basename.CDS");

			  	  &Rpred($Basename);
			  	  &SeqExtaction_Trans($Basename,$Tempfile);   
				 }

	else { print "Please Check the use option for RGPred Processing"; exit(); }

my $end_run = time();
my $run_time = $end_run - our $start_run;
print "\n\n!!RGPred running time for dataset $Infile is $run_time seconds!!\n";




		  ######################################## ******* SUBROUTINE  SECTION ******* ###########################################
				         ############################################################################

	sub Transdecoder()
	{
 				my $file = $Basename; my $basename = ''; my $index = ''; $file=~ s/\s//g; my $Tempfile1=$Tempfile;
			# CDS nucleotide and protein sequence extraction through Transdecoder
				system("TransDecoder.LongOrfs -t $file >/dev/null");
				system("TransDecoder.Predict -t $file >/dev/null");

			#Deletion of un-necessary files
				if($Tempfile1 eq "NO" ) { system("rm -rf *checkpoints *longorfs *.transdecoder_dir"); system("rm -f *.bed *.gff3 *cmds");  }	
				print "\n\n ########### 1. SUBROUTINE :Transdecoder CDS extraction: COMPLETED #################\n\n";sleep(2);
				
	}

	sub Augustus()
	{
				my ($file,$basename,$species) = @_; 
				system("augustus --species=$species --gff3=on --uniqueGeneId=true --codingseq=on $file >$basename.gff");
				system("./Exe/getAnnoFasta.pl $basename.gff");
		
	}

	sub Rpred()
		{
			 my $line='';my @arr_model=''; my @arr_model1=''; my @arr_model2=''; my @arr_model3=''; my @arr_model4=''; my $r = ''; my $basename=$Basename;
			 my @arr_model5='';  my $cutoff=$Score_cutoff; my $Pepfile=$basename.".PEP"; my $Tempfile1=$Tempfile; print "$Pepfile\n";
			# 1. sequence composition calculation and Vector creation 
				system("./Exe/aa_calc_0 $Pepfile $Pepfile.aafreq_out");
				system("./Exe/dipep_calc $Pepfile $Pepfile.dipep_out");
				system("./Exe/tripep_calc $Pepfile $Pepfile.tripep_out");
				system("./Exe/multipep_calc $Pepfile $Pepfile.mult_out");
				system("./Exe/hydro_calc $Pepfile $Pepfile.hdr_out");
				system("./Exe/cc_calc $Pepfile $Pepfile.cc_out");

				open (FA,"<$Pepfile.aafreq_out") or die "Unable to open $!"; my @arra=<FA>;
				open (FD,"<$Pepfile.dipep_out") or die "Unable to open $!"; my @arrd=<FD>;
				open (FT,"<$Pepfile.tripep_out") or die "Unable to open $!"; my @arrt=<FT>;
				open (FM,"<$Pepfile.mult_out") or die "Unable to open $!"; my @array=<FM>;
				open (FH,"<$Pepfile.hdr_out") or die "Unable to open $!"; my @arrh=<FH>;
				open (FC,"<$Pepfile.cc_out") or die "Unable to open $!"; my @arrcc=<FC>;
				open (FF, ">$Pepfile.all.dat") or die "Unable to open $!";
				  $r=0;
					foreach $line(@arra)
					{
						chomp $line; print FF "$line"; chomp $arrd[$r];print FF "$arrd[$r]";
						chomp $arrt[$r];print FF "$arrt[$r]"; chomp $array[$r];print FF "$array[$r]";
						chomp $arrh[$r];print FF "$arrh[$r]"; print FF "$arrcc[$r]";
						$r++;
					}
				close FD;close FM;close FH;close FC;close FT; close FA;close FF;
				print "\n\n ########### 3.SUBROUTINE :Physiochemical property calculation: COMPLETED #################\n\n"; sleep(2);

	
	      	 	# 2. SVM Classification
				system("./Exe/svm_classify $Pepfile.all.dat ./Exe/CNL_mixPara6_nbs_model471 CNL_Model_result");
				system("./Exe/svm_classify $Pepfile.all.dat ./Exe/TNL_mixPara6_nbs_model450 TNL_Model_result");
				system("./Exe/svm_classify $Pepfile.all.dat ./Exe/RLK_mixPara6_nbs_model450 RLK_Model_result");
				system("./Exe/svm_classify $Pepfile.all.dat ./Exe/RLP_mixPara6_nbs_model470 RLP_Model_result");
								
				open(FM0,"<CNL_Model_result") or die "failed to open file because of $!\n";
				open(FM1,"<TNL_Model_result") or die "failed to open file because of $!\n";
				open(FM2,"<RLK_Model_result") or die "failed to open file because of $!\n";
				open(FM3,"<RLP_Model_result") or die "failed to open file because of $!\n";
				
				@arr_model=<FM0>; @arr_model1=<FM1>; @arr_model2=<FM2>; @arr_model3=<FM3>;
				close FM0; close FM1;close FM2;close FM3;

			# 3. SVM Output Writing
			        `cat $Pepfile|grep "\>"|cut -d " " -f 1|sed 's/\>//' >$Pepfile.header`;
				open(FI,"<$Pepfile.header") or die "failed to open because $!\n"; @head=<FI>; close FI;
				open(OUT,">$Basename.RPredResult") or die "failed to open file because of $!\n";
				my $d='';my $h=0;my $e=0;  my $score="";
					foreach $d(@head)
					{ chomp $d; $score=0;
						if(($arr_model[$h]>$arr_model1[$h]) && ($arr_model[$h]> $arr_model2[$h]) && ($arr_model[$h]> $arr_model3[$h])&& $arr_model[$h]>=$cutoff)
							{ $score=$arr_model[$h]; chomp $score; print OUT "$d\t$score\tCNL\n";print "$d\t$score\tCNL\n";} 

						elsif(($arr_model1[$h]> $arr_model[$h]) && ($arr_model1[$h]> $arr_model2[$h]) && ($arr_model1[$h]> $arr_model3[$h])&& $arr_model1[$h]>=$cutoff)
							{ $score=$arr_model1[$h]; chomp $score; print OUT "$d\t$score\tTNL\n";print "$d\t$score\tTNL\n";} 

						elsif(($arr_model2[$h]> $arr_model[$h]) && ($arr_model2[$h]> $arr_model1[$h]) && ($arr_model2[$h]> $arr_model3[$h])&& $arr_model2[$h]>=$cutoff)
							{ $score=$arr_model2[$h]; chomp $score; print OUT "$d\t$score\tRLK\n";print "$d\t$score\tRLK\n";} 

						elsif(($arr_model3[$h]> $arr_model[$h]) && ($arr_model3[$h]> $arr_model1[$h]) && ($arr_model3[$h]> $arr_model2[$h])&& $arr_model3[$h]>=$cutoff)
							{ $score=$arr_model3[$h]; chomp $score; print OUT "$d\t$score\tRLP\n"; print "$d\t$score\tRLP\n";}

						else {print OUT "$d\t$score\tNonRGene\n";print "$d\t$score\tNonRGene\n";}
						$h++;
					}
				#Deletion of un-necessary files

				if($Tempfile eq "NO" ) { system("rm *_out" ); system("rm *Model_result"); system("rm *all.dat"); system("rm $Pepfile.header"); }
				print "\n\n ########### 4.SUBROUTINE :Rgene prediction: COMPLETED #################\n\n"; sleep(2);

		}

	sub SeqExtaction_Trans
		     {		
			#4. Extraction of Resistance gene nucleotide sequences
				 my $Cdsfile=$Basename.".CDS";  my $Pepfile=$Basename.".PEP"; my $Tempfile1=$Tempfile; my $Basename=$Basename;

				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "RLK"|cut -f 1|sort|uniq >$Cdsfile.RLK_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "RLP"|cut -f 1|sort|uniq >$Cdsfile.RLP_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "CNL"|cut -f 1|sort|uniq >$Cdsfile.CNL_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "TNL"|cut -f 1|sort|uniq >$Cdsfile.TNL_ids`;

				#`cat $Pepfile.RPredResult|grep -v -w "NonRGene"|cut -f 1|sort|uniq >$Pepfile.Selected_RGene_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|sort|uniq >$Cdsfile.RPredResult`; 


				`cat $Cdsfile|sed 's/^>/#>/'|sed 's/\$/%/'|sed 's/!^>\$//'|tr -d "\n"|tr "#" "\n"|sed 's/%/@/'|sed 's/%//g' >$Cdsfile"_singleline.fasta"`;

				`cat $Cdsfile"_singleline.fasta"|grep -f $Cdsfile.RLK_ids|sed 's/\|/\|RLK\|/'|tr "@" "\n" >$Cdsfile"_RLK_Ext.fn"`;
				`cat $Cdsfile"_singleline.fasta"|grep -f $Cdsfile.RLP_ids|sed 's/\|/\|RLP\|/'|tr "@" "\n" >$Cdsfile"_RLP_Ext.fn"`;
				`cat $Cdsfile"_singleline.fasta"|grep -f $Cdsfile.CNL_ids|sed 's/\|/\|CNL\|/'|tr "@" "\n" >$Cdsfile"_CNL_Ext.fn"`;
				`cat $Cdsfile"_singleline.fasta"|grep -f $Cdsfile.TNL_ids|sed 's/\|/\|TNL\|/'|tr "@" "\n" >$Cdsfile"_TNL_Ext.fn"`;

				`cat *_Ext.fn >$Cdsfile"_Ext_RGeneSeq.fn"`; `rm *_Ext.fn`;
				#`rm $Cdsfile"_CNL_Ext.fa $Cdsfile"_TNL_Ext.fa $Cdsfile"_RLK_Ext.fa $Cdsfile"_RLP_Ext.fa`;
				#`cat $Cdsfile"_singleline.fasta"|grep -f $Pepfile.Selected_RGene_ids|tr "@" "\n" >$Cdsfile"_Ext_RGeneSeq.fn"`;

				&SeqExtaction_Prot($Basename,$Tempfile); 


				#`cat $Pepfile|sed 's/^>/#>/'|sed 's/\$/%/'|sed 's/!^>\$//'|tr -d "\n"|tr "#" "\n"|sed 's/%/@/'|sed 's/%//g' >$Pepfile"_singleline.fasta"`;
				#`cat $Pepfile"_singleline.fasta"|grep -f $Pepfile.Selected_RGene_ids|tr "@" "\n" >$Pepfile"_Ext_RGeneSeq.fa"`;

				#Deletion of un-necessary files
				if($Tempfile eq "NO" ) { system("*transdecoder*"); }

				print "\n\n ########### 5.SUBROUTINE :Extraction of Rgene sequence from transcriptome : COMPLETED #################\n\n"; sleep(2);
				print "\n\n ###########!!!!!!! ******* RGPred pipeline successfully:COMPLETED *******!!!!!!! #################\n\n"; sleep(1);
		     }

	sub SeqExtaction_Prot
		     {		
			#4. Extraction of Resistance gene protein sequences
				 my $Pepfile=$Basename.".PEP"; my $Tempfile1=$Tempfile; $Basename=$Basename;

				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "RLK"|cut -f 1|sort|uniq >$Pepfile.RLK_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "RLP"|cut -f 1|sort|uniq >$Pepfile.RLP_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "CNL"|cut -f 1|sort|uniq >$Pepfile.CNL_ids`;
				`cat $Basename.RPredResult|grep -v -w "NonRGene"|grep "TNL"|cut -f 1|sort|uniq >$Pepfile.TNL_ids`;
				#`cat $Pepfile.RPredResult|grep -v -w "NonRGene"|cut -f 1|sort|uniq >$Pepfile.Selected_RGene_ids`;

				`cat $Basename.RPredResult|grep -v -w "NonRGene"|sort|uniq >$Pepfile.RPredResult`; 

				`cat $Pepfile|sed 's/^>/#>/'|sed 's/\$/%/'|sed 's/!^>\$//'|tr -d "\n"|tr "#" "\n"|sed 's/%/@/'|sed 's/%//g' >$Pepfile"_singleline.fasta"`;

				
				`cat $Pepfile"_singleline.fasta"|grep -f $Pepfile.RLK_ids|sed 's/\|/\|RLK\|/'|tr "@" "\n" >$Pepfile"_RLK_Ext.fa"`;
				`cat $Pepfile"_singleline.fasta"|grep -f $Pepfile.RLP_ids|sed 's/\|/\|RLP\|/'|tr "@" "\n" >$Pepfile"_RLP_Ext.fa"`;
				`cat $Pepfile"_singleline.fasta"|grep -f $Pepfile.CNL_ids|sed 's/\|/\|CNL\|/'|tr "@" "\n" >$Pepfile"_CNL_Ext.fa"`;
				`cat $Pepfile"_singleline.fasta"|grep -f $Pepfile.TNL_ids|sed 's/\|/\|TNL\|/'|tr "@" "\n" >$Pepfile"_TNL_Ext.fa"`;
				
				`cat *_Ext.fa >$Pepfile"_Ext_RGeneSeq.fa"`; `rm *_Ext.fa`;


				#Deletion of un-necessary files
				if($Tempfile eq "NO" ) { system("rm *_ids $Basename"); system("rm *.clean *.NCBI *singleline.fasta $Basename.RPredResult"); }

				print "\n\n ########### 5.SUBROUTINE :Extraction of Rgene sequence from proteome : COMPLETED #################\n\n"; sleep(2);
				print "\n\n ###########!!!!!!! ******* RGPred pipeline successfully:COMPLETED *******!!!!!!! #################\n\n"; sleep(1);

				system("/data/bioinfo/interproscan-5.27-66.0/interproscan.sh -i $Pepfile'_Ext_RGeneSeq.fa' -t p -dp -pa --goterms --iprlookup --cpu 5 -o $Pepfile'_Ext_RGeneSeq.iprscan'")







				
		     }

	sub Format_Ambiguity_removal
		     	{
				my ($file,$basename,$type) = @_;  my $index = ''; $_[0] = 1;
				open (IN, $file);
				open (OUT, ">$file.NCBI");	$index = 1;
					while (<IN>)
	 				{	s/\s+$//;
     	 					if ($_ =~ /^>/) { s/^>//; print OUT '>RPred',$index, '|', $_, "\n";	$index++; } else   { print OUT $_, "\n"; }
					} 
			  		 close(IN);close(OUT);

				open(UPI,"<$file.NCBI") or die "failed to open because $!\n"; my @UPI=<UPI>; close UPI;
				open(UPO,">$file.Ncbi_ambiguity_clean") or die "failed to open file because of $!\n"; 
      				foreach my $UPLI_line(@UPI)
  					{	#chomp $d;
						if($UPLI_line=~/^>/) {	print UPO "$UPLI_line"; } else { $UPLI_line =~ s/B|J|O|X|Z|U//g;	print UPO "$UPLI_line";}	
					}     close UPO;	
		        }


	sub filter_by_len
		{
				my ($file, $basename) = @_; 
				print "filtering of input file by length : $file $basename in progress!!!!\n"; sleep (2);
				system("./Exe/filter_by_length -o $file".".clean $file -n 50 -x 1000000");
		}

	sub filter_duplicate
		{
				my ($file, $basename) = @_; 
				print "filtering of duplicates in input file:$file $basename in progress!!!!\n"; sleep (2);
				system("./Exe/filter_duplicates -o $basename".".unduplicated $file".".clean ");
				`mv $basename"."unduplicated $basename`;
		}
