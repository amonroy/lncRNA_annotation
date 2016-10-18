import subprocess
import sys

file = sys.argv[1]

with open(file, 'r') as f:
	while True:
		seq_a = open('/tmp/temporary_seq_a', 'w')

		seq_a_header = f.readline()
		seq_a_seq = f.readline()
		
		seq_a.write(seq_a_header)
		seq_a.write(seq_a_seq)
		
		seq_a.close()
		
		seq_b = open('/tmp/temporary_seq_b', 'w')
		seq_b_header = f.readline()
		seq_b_seq = f.readline()

		seq_b.write(seq_b_header)
		seq_b.write(seq_b_seq)
		
		seq_b.close()

		proc = subprocess.Popen('water -asequence /tmp/seq_temporary_seq_a -bsequence /tmp/temporary_seq_b -gapextend 0.5 -gapopen 10.0 -aformat3 markx3 -auto -stdout >> /netscr/anais/foobar', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = proc.communicate()
		if err:
			print "Error with: %s & %s " % (seq_a_header.rstrip(), seq_b_header.rstrip())


	proc = subprocess.Popen('rm /tmp/temporary_seq_a', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	bproc = subprocess.Popen('rm /tmp/temporary_seq_b', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#water -asequence FBtr0340236-mel.fasta -bsequence FBtr0340236-sec.fasta -gapextend 0.5 -gapopen 10.0 -outfile FBtr0340236.water.markx10 -aformat3 markx10
