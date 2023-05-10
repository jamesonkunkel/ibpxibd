/*
    The Chromosome object is used by the Individuals to store genetic data. Mutations and recombination will occur within the Chromosome object.
*/

import java.util.Random;

public class Chromosome {

    //counts of QTL and neutral loci in the simulated genomes
    final int QTL_COUNT = 5;
    final int NEUTRAL_COUNT = 5;

    //length of a diploid chromosome with both QTL and neutral loci
    final int GENOME_LENGTH = QTL_COUNT + NEUTRAL_COUNT;

    //random generator for Gaussian and uniform allelic values
    static Random generator = new Random();

    //integer arrays for complementary pairs of haploid genome
    int[] genome1 = new int[GENOME_LENGTH];
    int[] genome2 = new int[GENOME_LENGTH];

    //randomly QTL and neutral loci with allelic values, Gaussian if QTL, uniform 1 or 0, if neutral, takes a standard deviation for the width of Gaussian kernel
    public static void fillGenome(int[] genome){
        for(int i = 0; i < genome.length; i++){
            int ran = generator.nextInt(2);

            if(ran == 0){
                genome[i] = 1;
            }else{
                genome[i] = -1;
            }
        }
    }

    public static void fillGenomeEmpty(int[] genome){
        for(int i = 0; i < genome.length; i++){
            int ran = generator.nextInt(2);

            if(ran == 0){
                genome[i] = 0;
            }else{
                genome[i] = 0;
            }
        }
    }

    //inserts a mutation into a specific loci for a specific genome
    public void insertMutation(int genome, int indexOfInsertion, int sd){

        int value = (int) (generator.nextGaussian() * sd);
        
        if(genome == 0){
            this.genome1[indexOfInsertion] = value;
        }else{
            this.genome2[indexOfInsertion] = value;
        }
    }
    
    //default contructor preloads both genomes with random integer values
    public Chromosome() {
        fillGenomeEmpty(this.genome1);
        fillGenomeEmpty(this.genome2);
    }

    //overloaded constructor that takes in two genomes as params
    public Chromosome(int[] genome1, int[] genome2){
        this.genome1 = genome1;
        this.genome2 = genome2; 
    }
}