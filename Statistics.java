/*
    Contains all the statistical output functions
*/

import java.util.ArrayList;

public class Statistics {

    //mean population phenotypic value, sum of all allelic values across all loci
    public void calculateMeanPhenotype(ArrayList<Individual> individuals){
        int total = 0;

        for(int i = 0; i < individuals.size(); i++){
            total += individuals.get(i).getPhenotype();
        }

        System.out.println("Mean phenotype is: " + total/individuals.size());
    }

    //calculates the phenotypic variance of the adult population
    public double calculatePhenotypicVariance(ArrayList<Individual> individuals){
        double total = 0;

        for(int i = 0; i < individuals.size(); i++){
            total += individuals.get(i).getPhenotype();
        }

        double mean = total / individuals.size();

        double sum = 0;
        for(int j = 0; j < individuals.size(); j++){
            sum += Math.pow((individuals.get(j).getPhenotype() - mean), 2);
        }

        double variance = sum/individuals.size();

        return variance;
    }

    //calculates the phenotypic variance of the adult population
    public double calculateGeneticVariance(ArrayList<Individual> individuals){
        double total = 0;

        for(int i = 0; i < individuals.size(); i++){
            total += individuals.get(i).getGeneticValue();
        }

        double mean = total / individuals.size();

        double sum = 0;
        for(int j = 0; j < individuals.size(); j++){
            sum += Math.pow((individuals.get(j).getGeneticValue() - mean), 2);
        }

        double variance = sum/individuals.size();

        return variance;
    }

    //default constructor method
    public Statistics() {

    }
}
