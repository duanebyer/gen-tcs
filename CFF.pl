#!/apps/bin/perl

#require "ctime.pl";

#open cff file
$filename = "CFF_DD_March2020.dat";
open (FILE, ">/vol0/pierre/Bureau/TCSGen/ChangementParam/$filename");


#Boucle en q2        
$num1 = 1;
$num2 = 10;

$PI = 3.14159265359;

for ($i = 0; $i <= 16; $i++) {

	$q2=(1.)+(0.5)*($i);

#Boucle en t        
	$num1t = 1;
	$num2t = 3;

	for ($j = 0; $j <= 16; $j++) {

		$t= -1*((0.)+(0.05)*($j));
		$t1= ((0.)+(0.05)*($j));

#Boucle en xi        
		$num1x = 1;
		$num2x = 2;

		for ($k = 0; $k <= 8; $k++) {

			
			$xi= (0.0526)+(0.08)*($k);
			
			$filenameCMD = "input.txt";
   			open (FILEcmd, ">/vol0/pierre/Bureau/TCSGen/ChangementParam/$filenameCMD");
    			print FILEcmd "3\n";
    			print FILEcmd "12\n";
    			print FILEcmd "34\n";
    			print FILEcmd "15\n";
    			print FILEcmd "1\n";
    			print FILEcmd "${q2}\n";
    			print FILEcmd "${t1}\n";
    			print FILEcmd "2\n";
    			print FILEcmd "1\n";
    			print FILEcmd "${xi}\n";
    			close (FILEcmd);
    			$cmd = "./dvcs < input.txt";
    			system($cmd);


			my $filename1 = 'cffs.dat';
			open(my $fh, '<:encoding(UTF-8)', $filename1)
			  or die "Could not open file '$filename1' $!";
 
			my $row = <$fh>;
			chomp $row;
			
			$dterm=0.0;
			
			my @vals = split(' ',$row);
			print FILE "${q2} ${t} ${vals[0]} ${vals[6]} ${vals[2]} ${vals[7]} ${vals[3]} ${vals[8]} ${vals[4]} ${dterm}\n";
			


		}
	}
}

close (FILE);
