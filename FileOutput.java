import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class FileOutput {

    public void outputVariances(double[] variances, String path){
        try {
            FileWriter writer = new FileWriter(path + "/phenVars.csv");

            writer.write("Generation,Phenotypic Variance\n");

            for(int i = 0; i < variances.length; i++){
                if(i == 0 | (i + 1) % 10 == 0){
                    writer.append(Integer.toString(i + 1) + "," + Double.toString(variances[i]));
                    writer.append("\n");
                }
            }

            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println(e.getLocalizedMessage());
        }
    }

    public void outputGenVariances(double[] variances, String path){
        try {
            FileWriter writer = new FileWriter(path + "/genVars.csv");

            writer.write("Generation,Genetic Variance\n");

            for(int i = 0; i < variances.length; i++){
                if(i == 0 | (i + 1) % 10 == 0){
                    writer.append(Integer.toString(i + 1) + "," + Double.toString(variances[i]));
                    writer.append("\n");
                }
            }

            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println(e.getLocalizedMessage());
        }
    }

    public void outputPhenotypeLocation(ArrayList<Individual> individuals, String path, int generation){
        try {
            File file = new File(path + "/phenLocation");
            file.mkdirs();
            FileWriter writer = new FileWriter(path + "/phenLocation/phenLocation" + Integer.toString(generation) + ".csv");

            writer.write("Individual,Phenotype,Genetic,GeneticNeut,X,Y,Offspring,Mat,Pat\n");

            for(int i = 0; i < individuals.size(); i++){
                Individual ind = individuals.get(i);

                writer.append(Integer.toString(i) + "," + Integer.toString(ind.getPhenotype()) + "," + Integer.toString(ind.getGeneticValue()) + "," + Integer.toString(ind.getGeneticValueNeut()) + "," + Double.toString(ind.getXCoord()) + "," + Double.toString(ind.getYCoord()) + "," + Integer.toString(ind.getOffspringCount()) + "," + Integer.toString(ind.getOffspringCountMat()) + "," + Integer.toString(ind.getOffspringCountPat()));
                writer.append("\n");
            }

            writer.flush();
            writer.close();
        } catch (Exception e) {
            System.out.println(e.getLocalizedMessage());
        }
    }

    public void outputQTLGenomes(ArrayList<Individual> individuals, String path, int generation){
        try {
            File file = new File(path + "/QTLgenomes");
            file.mkdirs();
            FileWriter writer = new FileWriter(path + "/QTLgenomes/genomes" + Integer.toString(generation) + ".csv");

            //writing headers for genome CSV file, Q = QTL, N = neutral locus
            for(int i = 0; i < 10; i++){
                writer.append("Q" + Integer.toString(i + 1) + ",");
            }
            writer.append("\n");


            for(int j = 0; j < individuals.size(); j++){
                for(int k = 0; k < 2; k++){
                    for(int l = 0; l < 5; l++){
                        if(k == 0){
                            writer.append(Integer.toString(individuals.get(j).getGenome1()[l]) + ",");
                        }else{
                            writer.append(Integer.toString(individuals.get(j).getGenome2()[l]) + ",");
                        }
                    }
                }
                writer.append("\n");
            }

            writer.flush();
            writer.close();
        } catch (Exception e) {
            System.out.println(e.getLocalizedMessage());
        }
    }

    public void outputNeutralGenomes(ArrayList<Individual> individuals, String path, int generation){
        try {
            File file = new File(path + "/neutGenomes");
            file.mkdirs();
            FileWriter writer = new FileWriter(path + "/neutGenomes/genomes" + Integer.toString(generation) + ".csv");

            //writing headers for genome CSV file, Q = QTL, N = neutral locus
            for(int i = 0; i < 10; i++){
                writer.append("Q" + Integer.toString(i + 1) + ",");
            }
            writer.append("\n");

            for(int j = 0; j < individuals.size(); j++){
                for(int k = 0; k < 2; k++){
                    for(int l = 5; l < 10; l++){
                        if(k == 0){
                            writer.append(Integer.toString(individuals.get(j).getGenome1()[l]) + ",");
                        }else{
                            writer.append(Integer.toString(individuals.get(j).getGenome2()[l]) + ",");
                        }
                    }
                }
                writer.append("\n");
            }

            writer.flush();
            writer.close();
        } catch (Exception e) {
            System.out.println(e.getLocalizedMessage());
        }
    }

    public FileOutput(){

    }
}
