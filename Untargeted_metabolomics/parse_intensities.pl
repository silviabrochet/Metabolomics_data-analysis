#! /usr/bin/perl -w

##Example usage:
# perl parse_intensities.pl bee_May2017_test_strain_specificity.txt > my_parsed_filt.txt

#Sub-routine for trimming
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

#Read data from both sheets
open (FILE1, $ARGV[0]) or
  die "Can't open $ARGV[0]: $!";
my @file1 = <FILE1>;
my $header1 = shift @file1;
my @split_header1 = split("\t",$header1);
my @sample_ids1 = splice @split_header1,3;

close FILE1;

open (FILE2, $ARGV[1]) or
    die "Can't open $ARGV[1]:$!";
my @file2 = <FILE2>;
my $header2 = shift @file2;
my @split_header2 = split("\t",$header2);
my @sample_ids2 = splice @split_header2,3;
close FILE2;

my @all_sample_ids = (@sample_ids1,@sample_ids2);
my %tech_id;
foreach (@all_sample_ids) {
    chomp;
    my @split_id = split("\/",$_);
    my $trim_id = sprintf("%.0f",trim($split_id[0]));
    my $trim_name = trim($split_id[1]);
    $tech_id{$trim_id} = $trim_name;
}
my $nb_tech_id = keys %tech_id;

#Associate each sample with the two technical replicate ids

my %sample_id;
foreach (keys %tech_id){
    my $sample = $tech_id{$_};
    $sample_id{$sample} = [] unless exists $sample_id{$sample};
    push @{$sample_id{$sample}},$_;
}
my @sorted_samples = sort keys %sample_id;
my $nb_samples = @sorted_samples;

#Associate each ion-id with an array of ion intensities, corresponding to each of the samples

my $i = 0;
my $nb_ions = @file1;
my %ion_line;
while($i<$nb_ions) {
    my @split_file1 = split("\t",$file1[$i]);
    my @measurements1 = splice @split_file1,3;
    my @split_file2 = split("\t",$file2[$i]);
    my @measurements2 = splice @split_file2,3;
    my @all_measurements = (@measurements1,@measurements2);
    my $nb = @all_measurements;
    my $ion_id = $split_file1[0];
    $ion_line{$ion_id} = \@all_measurements;
    ++$i;
}
my @sorted_ions = sort({$a <=> $b} keys %ion_line);

#Next, t-transform the data, so that each sample-id is associated with a list of ion-intensities for all detected ions. 

my %sample_measurements ;
my $j = 0;
while ($j < $nb_tech_id) {
    $sample_measurements {$j} = [];
    foreach my $ion(@sorted_ions) {
	push @{$sample_measurements {$j}},$ion_line{$ion}[$j];
    }
    ++$j;
}

#Calculate the mean-values of the technical replicates for each sample

my %sample_measurements_mean;
foreach my $sample(@sorted_samples) {
    my @tech_ids = sort @{$sample_id{$sample}};
    my $k = 0;  
    my @sample_ion_means=();
    while ($k < $nb_ions) {
	my $mean_value = ($sample_measurements{($tech_ids[0]-1)}[$k]+$sample_measurements{($tech_ids[1]-1)}[$k])/2;
	push @sample_ion_means,$mean_value;
	++$k;
    }
    $sample_measurements_mean{$sample} = \@sample_ion_means;
}

#Print an output file containing the ordered sample names in column 1, the sample group as a factor in column 2, and the mean sample measurements for each ion in the remaining columns. Sample names were not very convenient for inferring group-names, the script will need to be adjusted for future experiments for that block of code.

#print "Sample\tGroup\t",join("\t",@sorted_ions),"\n";
my %sample_groups;
foreach(@sorted_samples) {
    my @split_name = split("_",$_);
    my $sample_group;
    if ($split_name[0] eq "F5") {	
	$sample_group = join("_",@split_name[0..3]);
    }
    elsif ($split_name[0] eq "Neg") {
	$sample_group = join("_",@split_name[0..2]);
    }
    elsif ($split_name[0] eq "H2O") {
	$sample_group = $split_name[0]."_".$split_name[1];
    }
    elsif ($split_name[0] eq "cfMRS") {
	$sample_group = $split_name[0];
    }
    else {
	$sample_group = $_;
    }
    $sample_groups{$sample_group} =1;
#    print $_,"\t",$sample_group,"\t",join("\t",@{$sample_measurements_mean{$_}}),"\n";
}


foreach (keys %sample_groups) {
    if ($_ =~ /F5/ && $_ =~ /T16/) {
	print $_,"\n";
    }
}
