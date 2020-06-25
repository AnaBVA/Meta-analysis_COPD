#!/usr/bin/perl
use strict;
use Spreadsheet::WriteExcel;

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ DESCRIPTION /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\ Script that parses a flat file containing GEO data. This type of files clusters samples as follows:
#/\/\/\     No. Description
#/\/\/\	    Extended description (optional)
#/\/\/\     Organism
#/\/\/\     Type/Source name/number of:DataSets,Series,Related-Platforms,Samples (optional)
#/\/\/\     Platform Series No_samples(optional)
#/\/\/\     FTP GEO/SRA web-link
#/\/\/\     Dataset/Series/Sample 	Accession	 ID
#/\/\/\ Input: file containing GEO data (e.g. copd.txt)
#/\/\/\ Output: file containing GEO data as: 1) table format [tab delimited] (e.g. copd_tab.txt) and 
#/\/\/\         2) excel file (e.g. copd.xls)
#/\/\/\ Script name: parsecopd.pl 
#/\/\/\ Execution: perl parsecopd.pl infilename
#/\/\/\ *** This script requires a perl module. More information at bottom ***
#/\/\/\ Written by: Alejandra Zayas <ezayas@lcg.unam.mx> & Orlando Santillan <osantillan@lcg.unam.mx>. August, 2016
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

#/\/\/\/\/\/\/\/\/\/\/\/\/\ VARIABLE DECLARATION /\/\/\/\/\/\/\/\/\/\/\/\/\
my ($count, $report, @unique, $num, $desc, $ext_desc, $species, $num_dataset, $num_series, $num_relatedplatform, $num_samples);
my ($typesource, $platform, $series, $dataset, $geo, $ftp, $type, $accession, $id, $info);

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ MAIN SCRIPT /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
if(scalar @ARGV != 0){
 if(-e $ARGV[0]){
  &countentries($ARGV[0]);
  &parsefile($ARGV[0]);
  &makexcel($ARGV[0]);
 }else{
   print "  Given input file does NOT exist: ", $ARGV[0], "\n";
 }#end if-else. Verifies infile exists
}else{
 print "  Missing argument: input file name \n  Run: perl parsecopd.pl infile \n";
}#end if-else. Verifies infile was given
exit;

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ SUBROUTINES /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
sub countentries{
 open INFILE, "<", $ARGV[0]||die "  Can't read ", $ARGV[0], "\n ", $!, "\n";
 while(<INFILE>){
  chomp;
  if($_ =~ /^(\d+)\..+/){
   $count++;
   $report = $1;
  }elsif($_ =~ /.+\tID: (\d+)/){
   &unique($1);
  }#end if
 }#end while. Loops through INFILE
 close INFILE;
 if($count == $report){
  print "  ", $count, " entries found \n";
 }else{
  print "  ", $count, " found, but ", $report, "entries reported \n";
 }#end if-else. Prints number of entries found
 undef $report;
 if(scalar @unique == $count){
  print "  ", scalar @unique, " unique IDs found \n";
 }else{
  print "  ", scalar @unique, " unique IDs found, but", $count, "entries on input file \n";
 }#end if-else. Prints number of unique IDs found
 return $count;
}#end subroutine countentries

sub unique{
 my $element; my $trick = 0;
 foreach $element (@unique){
  if($element == $1){
   $trick++;
  }#end if
 }#end foreach
 if($trick == 0){
  push @unique, $1;
 }#end if
 return @unique;
}#end subroutine unique

sub parsefile{
 $ARGV[0] =~ /(.+)\.(.+)/;
 open OUTFILE, ">", $1."_tab.".$2 || die "  Can't write ", $1."_tab.".$2, "\n ", $!, "\n";
 print OUTFILE join "\t", 'No','Description','Extended Description','Organism','Type/Source','No. Datasets','No. Series','No. Rel. Platforms','No. Samples','Platform','Series','Dataset','GEO/SRA','FTP','Type','Accession',"ID\n";
 open INFILE, "<", $ARGV[0] || die "  Can't read ", $ARGV[0], "\n ", $!, "\n";
 $num = $desc = $ext_desc = $species = $num_dataset = $num_series = $num_relatedplatform = $num_samples = 'NA';
 $typesource = $platform = $series = $dataset = $geo = $ftp = $type = $accession = $id = 'NA'; $info = "";
 while(<INFILE>){
  chomp;
  if($_ =~ /^(\d+)\. (.+)/){
   $num = $1;
   $desc = $2;
   $info = join "\t", $num, $desc;
   
  }elsif($_ =~ /^(\w+.+|\(Submitter supplied\).+)\.$/ || $_ =~ /^(\(Submitter supplied\).+)/){
   $ext_desc = $1;
   
  }elsif($_ =~ /^Organism:\t(.+)/){
   $species = $1; 
   
  }elsif($_ =~ /^Type:|Source name:|\d+ DataSets|\d+ related Platforms/){
   if($_ =~ /^Type:\t{2}(.+)/){
    $typesource = $1;
   }elsif($_ =~ /^Source name:\t(.+)/){
    $typesource = $1;
   }elsif($_ =~ /^(\d+) DataSets (\d+) Series (\d+) Related Platforms (\d+) Samples$/){
    $num_dataset = $1;
    $num_series = $2;
    $num_relatedplatform = $3;
    $num_samples = $4;
   }elsif($_ =~ /^(\d+) related Platforms (\d+) Samples/){
    $num_relatedplatform = $1;
    $num_samples = $2;
   }#end if-elsif
   
  }elsif($_ =~ /^Platform/){
   if($_ =~ /^Platform: (.+) Series: (.+) (\d+) Samples$/){
    $platform = $1;
    $series = $2; 
    $num_samples = $3;
    $dataset = 'NA';
   }elsif($_ =~ /^Platform: (.+) (\d+) Samples$/){
    $platform = $1;
    $series = 'NA'; 
    $num_samples = $2;
    $dataset = 'NA';
   }elsif($_ =~ /^Platform: (.+) Series: (\w+\d+)$/){
    $platform = $1;
    $series = $2; 
    $num_samples = 'NA';
    $dataset = 'NA';
   }elsif($_ =~ /^Platform: (.+) Series: (.+) Dataset: (.+)$/){
    $platform = $1;
    $series = $2;
    $num_samples = 'NA';
    $dataset = $3;
   }elsif($_ =~ /^(Platform)\tAccession: (.+)\tID: (.+)/){
    $type = $1;
    $accession = $2;
    $id = $3;
   }#end if-elsif
   
  }elsif($_ =~ /^Dataset:/){
   if($_ =~ /^Dataset: (.+) Platform: (.+) (\d+) Samples$/){
    $platform = $2;
    $series = 'NA';
    $num_samples = $3;
    $dataset = $1;
   }#end if
   
  }elsif($_ =~ /^FTP download:/){
   if($_ =~ /^FTP download: GEO (ftp:\/\/ftp\..+)/){
    $geo = 'NA';
    $ftp = $1;
   }elsif($_ =~ /^FTP download: GEO \((.+)\) (.+)/){
    $geo = $1;
    $ftp = $2;
   }elsif($_ =~ /^FTP download: SRA (.+) (.+)/){
    $geo = $1;
    $ftp = $2;
   }#end if-elsif
   
  }elsif($_ =~ /^(DataSet|Series|Sample)\t{2}Accession: (.+)\tID: (.+)/){
   $type = $1;
   $accession = $2;
   $id = $3;
  }#end if-elsif
  
  &printoutfile($accession, $id);
 }#end while
 close INFILE;
 close OUTFILE;
}#end subroutine parsefile

sub printoutfile{
 if($accession ne 'NA' && $id ne 'NA'){
  $info = join "\t", $info, $ext_desc, $species, $typesource, $num_dataset, $num_series, $num_relatedplatform, $num_samples, $platform, $series, $dataset, $geo, $ftp, $type, $accession, $id;
  print OUTFILE $info, "\n";
  $info = ""; 
  $num = $desc = $ext_desc = $species = $num_dataset = $num_series = $num_relatedplatform = $num_samples = 'NA';
  $typesource = $platform = $series = $dataset = $geo = $ftp = $type = $accession = $id = 'NA';
 }#end if
}#end subroutine printoutfile

sub makexcel{
 my ($txtfile, $excelfile, $workbook, $worksheet, $formathead, $formatbody, $formatdescrip, $formatSp);
 my $row = -1; my (@line, $col);
 $ARGV[0] =~ /(.+)\.(.+)/;
 $txtfile = $1."_tab.".$2;
 $excelfile = $1.".xls";
 $workbook = Spreadsheet::WriteExcel -> new($excelfile);
 $worksheet = $workbook -> add_worksheet();
 $formathead = $workbook -> add_format();
 $formathead -> set_bold(); $formathead -> set_color('blue'); $formathead -> set_align('center');
 $formatbody = $workbook -> add_format();
 $formatbody -> set_align('center');
 $formatdescrip = $workbook -> add_format();
 $formatdescrip -> set_align('left');
 $formatSp = $workbook -> add_format();
 $formatSp -> set_italic(); $formatSp -> set_align('center');
 open TXT, "<", $txtfile || die " Can't read ", $txtfile, "\n ", $!, "\n";
 while(<TXT>){
  chomp;
  @line = split /\t/;
  $row++;
  for($col = 0; $col <= $#line; $col++){
   if($row == 0){
    $worksheet -> write($row, $col, $line[$col], $formathead);
   }else{
    if($col == 1 || $col == 2 || $col == 4 || $col == 14){
     $worksheet -> write($row, $col, $line[$col], $formatdescrip);
    }elsif($col == 3){
     $worksheet -> write($row, $col, $line[$col], $formatSp);
    }elsif($col == $#line){
     $worksheet -> write_number($row, $col, $line[$col], $formatbody);
    }else{
     $worksheet -> write($row, $col, $line[$col], $formatbody);
    }#end if-elsif-else
   }#end if-else
  }#end for
 }#end while. Loops through TXT file
 close TXT;
}#end subroutine makexcel

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ ***PERL MODULES CPAN*** /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#/\/\/\ Before running the script for the first time, type into a terminal:
#/\/\/\ 	sudo cpan App::cpanminus
#/\/\/\		sudo cpanm Spreadsheet::WriteExcel
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
