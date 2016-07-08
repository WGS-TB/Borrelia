library(Rsubread)
buildindex(basename ="ref1_index",reference ="/home/ali/Desktop/reference/ref2")
align(index="ref1_index",readfile1 = "/home/ali/Desktop/new.simulated data/seperated/AB.7.3_R1.fq 2",readfile2 =  "/home/ali/Desktop/new.simulated data/seperated/AB.7.3_R2.fq 2",output_file = "/home/ali/Desktop/new.simulated data/seperated/A.bam")