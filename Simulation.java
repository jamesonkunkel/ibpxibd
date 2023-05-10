import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.SwingUtilities;
import javax.swing.JFrame;

/*
    Description: Simulation objects represent individual replications of some simulation
    parameter set. They handle all the logic needed for handling the simulations. 
*/

public class Simulation {

    Random generator = new Random();
    Statistics stats = new Statistics();
    FileOutput output = new FileOutput();

    
    final static int MAX_DIM_X = 50;
    final static int MAX_DIM_Y = 50;

    //base path which depends on the scheme
    private String path;

    //random ID of simulation
    private int seed;

    //scheme of Simulation, 1=PAN, 2=IBP, 3=IBD, 4=IBPXIBD
    private int scheme;

    private int simLen;
    private int popSize;
    private double mutRate;

    private int genomeLength = 0;
    private int mutCount = 0;

    //s.d. of assortative mating function
    private double assortWidth;

    //s.d. of random Gaussian deviates for offspring dispersal
    private double disperseWidth;

    //s.d. of spatial mating function
    private double spatialWidth;

    //ArrayLists of currently mating Individual objects and offspring
    ArrayList<Individual> individuals = new ArrayList<Individual>();
    ArrayList<Individual> nextIndividuals = new ArrayList<Individual>();

    Landscape landscape = new Landscape(MAX_DIM_X, MAX_DIM_Y, individuals);

    //returns random poisson variable from distribution with shape lambda
    public int getPoisson(double lambda) {
        double L = Math.exp(-lambda);
        double p = 1.0;
        int k = 0;
    
        do {
        k++;
        p *= Math.random();
        } while (p > L);
    
        return k - 1;
    }

    //calcultes p.d.f for normal distrib.
    public double dnorm(double value, double mean, double sd){
        return (1/(sd*Math.sqrt(2*Math.PI))*Math.exp(-(Math.pow((value - mean), 2))/(2*Math.pow(sd, 2))));
    }

    //gets current date in dd/mm/yyyy format
    public String getDate(){
        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter format = DateTimeFormatter.ofPattern("MM-dd-yyyy");
        String formattedDate = now.format(format);

        return formattedDate;
    }

    //formats the output folder for a single simulation run
    public String formatOutputFolderName(int scheme, int seed, int popSize, double mutRate, double assortWidth, double disperseWidth, double spatialWidth){
        switch (scheme) {
            case 1:
                return getDate() + "_" + seed + "_popSize_" + popSize + "_mutRate_" + mutRate;
            case 2: 
                return getDate() + "_" + seed + "_popSize_" + popSize + "_mutRate_" + mutRate + "_assortWidth_" + assortWidth;
            case 3: 
                return getDate() + "_" + seed + "_popSize_" + popSize + "_mutRate_" + mutRate + "_disperseWidth_" + disperseWidth + "_spatialWidth_" + spatialWidth;
            case 4:
                return getDate() + "_" + seed + "_popSize_" + popSize + "_mutRate_" + mutRate + "_assortWidth_" + assortWidth + "_disperseWidth_" + disperseWidth + "_spatialWidth_" + spatialWidth;
            default:
                return "";
        }
    }

    //writes directory path that will store outputs from a simulation run
    public void makeOutputDirs(){
        switch (scheme) {
            case 1:
                this.path = this.path + "N=" + this.popSize + "/mutRate=" + this.mutRate + "/";
                break;
            case 2: 
                this.path = this.path + "N=" + this.popSize + "/mutRate=" + this.mutRate + "/assortWidth=" + this.assortWidth + "/";
                break;
            case 3: 
                this.path = this.path + "N=" + this.popSize + "/mutRate=" + this.mutRate + "/disperseWidth=" + this.disperseWidth + "/spatialWidth=" + this.spatialWidth + "/";
                break;
            case 4: 
                this.path = this.path + "N=" + this.popSize + "/mutRate=" + this.mutRate + "/assortWidth=" + this.assortWidth + "/disperseWidth=" + this.disperseWidth + "/spatialWidth=" + this.spatialWidth + "/";
                break;
            default:
                break;
        }

        this.path = this.path + formatOutputFolderName(this.scheme, this.seed, this.popSize, this.mutRate, this.assortWidth, this.disperseWidth, this.spatialWidth);

        File file = new File(this.path);
        file.mkdirs();
        System.out.println(this.path);
    }

    //generates n initial Individual objects
    public void addPopulation(int popSize){
        for(int i = 0; i < popSize; i++){
            individuals.add(new Individual(MAX_DIM_X, MAX_DIM_Y));
        }
    }

    //returns a 2 * GENOME_LENGTH matrix containing the freely recombined genome of a single individual
    public int[][] freelyRecombineGenome(int indexOfIndividual){
        //individual to perform free recombination on
        Individual individual = individuals.get(indexOfIndividual);

        //genomes of the individual
        int[] genome1 = individual.getGenome1();
        int[] genome2 = individual.getGenome2();

        // 2 x GENOME_LENGTH matrix that will store recombined genomes - row 1 is recombined genome 1, row 2 is recombined genome 2
        int[][] genomeMatrix = new int[2][individual.getGenomeLength()];

        //pick 0 or 1, if 0, do not swap alleles at given locus, otherwise recombine them (50:50 odds)
        for(int i = 0; i < individual.getGenomeLength(); i++){
            int ran = (int) (Math.random() * 2);

            if(ran == 0){
                genomeMatrix[0][i] = genome1[i];
                genomeMatrix[1][i] = genome2[i];
            }else{
                genomeMatrix[0][i] = genome2[i];
                genomeMatrix[1][i] = genome1[i];
            }
        }

    //return the recombined gametes in the form of an integer matrix
    return genomeMatrix;
    }

    //draw mutations and insert them into breeding population
    public void drawMutations(){
        if(this.mutCount == 0){
            this.genomeLength = (individuals.get(0).getGenomeLength());
            int mutDim = this.popSize * this.genomeLength * 2;
            this.mutCount = (int) (mutDim * this.mutRate);
        }

        int mutCountDeviate = getPoisson(1);
        int genMutCount = this.mutCount + mutCountDeviate;

        for(int i = 0; i < genMutCount; i++){
            int indexOfIndividual = (int) (Math.random() * this.popSize);
            int indexOfInsertion = (int) (Math.random() * this.genomeLength);
            int genomeOfIndividual = (int) (Math.random() * 2);
            
            individuals.get(indexOfIndividual).chromosome.insertMutation(genomeOfIndividual, indexOfInsertion, 3);
         }
    }

    //set location of offspring based on random deviate around maternal position
    public double[] disperseOffspring(Individual mother){
        double motherX = mother.getXCoord();
        double motherY = mother.getYCoord();
        
        double offspringX;
        double offspringY;

        do{
            offspringX = motherX + (generator.nextGaussian() * this.disperseWidth);
            offspringY = motherY + (generator.nextGaussian() * this.disperseWidth);
        }while(offspringX > landscape.getDimX() | offspringX < 0 | offspringY > landscape.getDimY() | offspringY < 0);

        double[] offspringCoords = {offspringX, offspringY};

        return offspringCoords;
    }

    public void changeOffspringCounts(int parent1Index, int parent2Index){
        Individual parent1 = individuals.get(parent1Index);
        Individual parent2 = individuals.get(parent2Index);

        parent1.incrementOffspringCount();
        parent2.incrementOffspringCount();
        parent1.incrementOffspringCountMat();
        parent2.incrementOffspringCountPat();
    }

    //mate two individuals to produce a new individual with free recombination
    public void mate(int parent1Index, int parent2Index){
        //increment offspring counts of both parents or just once if parent selfs
        changeOffspringCounts(parent1Index, parent2Index);

        //randomly select a genome from both individuals
        int genomeSelect1 = (int) (Math.random() * 2); 
        int genomeSelect2 = (int) (Math.random() * 2);

        //get the free recombination matrices of the two individuals
        int[][] parent1GenomeMatrix = freelyRecombineGenome(parent1Index);
        int[][] parent2GenomeMatrix = freelyRecombineGenome(parent2Index);

        //if genomeSelect is 0 then pass on genome 1, otherwise pass on genome 2
        int[] parent1Genome = genomeSelect1 == 0 ? parent1GenomeMatrix[0] : parent1GenomeMatrix[1];
        int[] parent2Genome = genomeSelect2 == 0 ? parent2GenomeMatrix[0] : parent2GenomeMatrix[1];
    
        double[] offSpringCoords = disperseOffspring(individuals.get(parent1Index));

        //add the offspring to the offspring ArrayList
        nextIndividuals.add(new Individual(parent1Genome, parent2Genome, offSpringCoords[0], offSpringCoords[1]));
    }

    //mates two individuals without recombination
    public void mateNoRecombination(int parent1Index, int parent2Index){
        //randomly select a genome from both individuals
        int genomeSelect1 = (int) (Math.random() * 2); 
        int genomeSelect2 = (int) (Math.random() * 2);

        //if genomeSelect is 0 then pass on genome 1, otherwise pass on genome 2
       int[] parent1Genome = genomeSelect1 == 0 ? individuals.get(parent1Index).getGenome1() : individuals.get(parent1Index).getGenome2();
       int[] parent2Genome = genomeSelect2 == 0 ? individuals.get(parent2Index).getGenome1() : individuals.get(parent2Index).getGenome2();

       double[] offSpringCoords = disperseOffspring(individuals.get(parent1Index));

       //add the offspring to the offspring ArrayList
       nextIndividuals.add(new Individual(parent1Genome, parent2Genome, offSpringCoords[0], offSpringCoords[1]));
   }

   //generates popSize offpsring for a single generation with random mating
   public void randomMating(int popSize){    
        for(int i = 0; i < this.popSize; i++){
            //randomly choose parents to mate POP_SIZE times 
            int parent1Index = (int) (Math.random() * this.popSize);
            int parent2Index = (int) (Math.random() * this.popSize);

            mate(parent1Index, parent2Index);
            //System.out.println("Ind " + parent1Index + " mates with ind " + parent2Index + " with value " + individuals.get(parent1Index).getPhenotype() + " and " + individuals.get(parent2Index).getPhenotype());
        }
    }

    //mates individuals assortatively based on weights from normal p.d.f
    public void assortativeMating(int popSize){
        for(int i = 0; i < popSize; i++){
            //randomly select mother
            int indexOfMother = (int) (Math.random() * popSize);

            //sum of weights which are pulled from normal p.d.f
            double weightSum = 0.0;

            for(int j = 0; j < popSize; j++){
                weightSum += dnorm(individuals.get(j).getPhenotype(), individuals.get(indexOfMother).getPhenotype(), this.assortWidth);
            }

            //index of father
            int indexOfFather = 0;

            //random value for weighted random draw of father
            double ran = Math.random() * weightSum;

            //sums weights for father selection
            double countWeight = 0.0;
            
            //sum weights until the index where the sum exceeds the random value then use that index as the father
            for(int k = 0; k < popSize; k++){
                countWeight += dnorm(individuals.get(k).getPhenotype(), individuals.get(indexOfMother).getPhenotype(), this.assortWidth);

                if(countWeight >= ran){
                    indexOfFather = k;
                    mate(indexOfMother, indexOfFather);
                    break;
                }
            }
            //System.out.println("Ind " + indexOfMother + " mates with ind " + indexOfFather + " with value " + individuals.get(indexOfMother).getPhenotype() + " and " + individuals.get(indexOfFather).getPhenotype());
        }
    }

    //generates POP_SIZE offpsring for a single generation with declining probability over distance between potential mates
    //when considering evaluation of exp(-(spatialWidth) * distance), larger spatialWidth values enforce stronger spatial mating
    public void spatialMating(int popSize, ArrayList<Individual> inds){
        double[] mateDistances = new double[popSize];
        
        for(int i = 0; i < popSize; i++){
            //randomly choose parents to mate POP_SIZE times 
             int indexOfMother = (int) (Math.random() * popSize);
            //int indexOfMother = i;

            double[] inverseDistances = new double[popSize];
            double weightSum = 0;

            for(int j = 0; j < popSize; j++){
                double mateDistance = landscape.getDistance(inds.get(indexOfMother), inds.get(j));
                double inverseDistance = Math.exp(-(this.spatialWidth)*mateDistance);

                inverseDistances[j] = inverseDistance;

                weightSum += inverseDistance;
            }

            //index of father
            int indexOfFather = 0;

            //random value for weighted random draw of father
            double ran = Math.random() * weightSum;

            //sums weights for father selection
            double countWeight = 0.0;

            for(int k = 0; k < popSize; k++){
                countWeight += inverseDistances[k];

                if(countWeight >= ran){
                    indexOfFather = k;
                    mate(indexOfMother, indexOfFather);
                    mateDistances[i] = landscape.getDistance(inds.get(indexOfMother), inds.get(indexOfFather));
                    break;
                }
            }
        }
    }

    //combines assortative mating with spatial mating pattern
    public void spatialAssortativeMating(){

        double[] mateDistances = new double[popSize];

        for(int i = 0; i < this.popSize; i++){
            //randomly choose parents to mate POP_SIZE times 
            // int indexOfMother = i;
            int indexOfMother = (int) (Math.random() * popSize);

            double[] productWeights = new double[this.popSize];
            double weightSum = 0;

            for(int j = 0; j < this.popSize; j++){
                double mateDistance = landscape.getDistance(this.individuals.get(indexOfMother), this.individuals.get(j));
                double inverseDistance = Math.exp(-(this.spatialWidth) * mateDistance);
                //System.out.println("Distance weight: " + inverseDistance);
                
                double assortativeWeight = dnorm(individuals.get(j).getPhenotype(), individuals.get(indexOfMother).getPhenotype(), this.assortWidth);
                //System.out.println("Assort weight: " + assortativeWeight);


                productWeights[j] = inverseDistance * assortativeWeight;
                //System.out.println("Product: " + productWeights[j]);

                weightSum += productWeights[j];
            }

            //index of father
            int indexOfFather = 0;

            //random value for weighted random draw of father
            double ran = Math.random() * weightSum;

            //sums weights for father selection
            double countWeight = 0;

            for(int k = 0; k < this.popSize; k++){
                countWeight += productWeights[k];

                if(countWeight >= ran){
                    indexOfFather = k;
                    mate(indexOfMother, indexOfFather);
                    mateDistances[i] = landscape.getDistance(individuals.get(indexOfMother), individuals.get(indexOfFather));
                    break;
                }
            }
        }
    }

    //transfers nextIndividuals into the individuals ArrayList in preparation for the next generation
    public void offspringToAdult(){
        //clear out current individuals after they have mated
        individuals.clear();
        
        //loop and add offspring into the adult ArrayList
        for(int i = 0; i < this.popSize; i++){
            individuals.add(nextIndividuals.get(i));
        }

        //clear out past offspring ArrayList
        nextIndividuals.clear();
    }

    //creates the visualizer frame and displays it
    private void createAndShowGUI() {
        System.out.println("Created GUI on EDT? "+
        SwingUtilities.isEventDispatchThread());

        JFrame f = new JFrame("Spatial Simulation Demo");
        f.setResizable(false);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.add(landscape);
        f.pack();
        f.setVisible(true);
    }

    //starts a simulation object
    public void runSimulationPAN() throws InterruptedException{
        //the main simulation loop will go here
        System.out.println("PAN simulation " + seed + " begins.");
        createAndShowGUI();
        addPopulation(popSize);

        int time = 0;
        int burnIn = 20000;

        double[] variances = new double[this.simLen];
        double[] genVars = new double[this.simLen];

        System.out.println("Burning in...");

        //burn-in loop
        while(time < burnIn){
            drawMutations();
            randomMating(popSize);
            offspringToAdult();

            time++;
        }

        System.out.println("Post burn-in...");

        //post burn-in
        while(time < (burnIn + simLen)){
            variances[time - burnIn] = stats.calculatePhenotypicVariance(individuals);
            genVars[time - burnIn] = stats.calculateGeneticVariance(individuals);

            drawMutations();
            randomMating(popSize);

            if(time % 10 == 0 | time == 0){
                output.outputPhenotypeLocation(individuals, path, time - burnIn);
                output.outputQTLGenomes(individuals, path, time - burnIn);
                output.outputNeutralGenomes(individuals, path, time - burnIn);
            }

            landscape.setIndividuals(nextIndividuals);
            offspringToAdult();

            time++;
        }

        output.outputVariances(variances, this.path);
        output.outputGenVariances(genVars, path);
    }

    //starts a simulation object
    public void runSimulationIBP() throws InterruptedException{
        //the main simulation loop will go here
        System.out.println("IBP simulation " + seed + " begins.");
        createAndShowGUI();
        addPopulation(popSize);

        int time = 0;
        int burnIn = 10000;
        double[] variances = new double[this.simLen];
        double[] genVars = new double[this.simLen];

        System.out.println("Burning in...");

        //burn-in loop
        while(time < burnIn){
            drawMutations();
            randomMating(popSize);
            offspringToAdult();

            time++;
        }

        System.out.println("Post burn-in...");

        while(time < (burnIn + simLen)){
            variances[time - burnIn] = stats.calculatePhenotypicVariance(individuals);
            genVars[time - burnIn] = stats.calculateGeneticVariance(individuals);
            
            drawMutations();
            assortativeMating(popSize);

            if(time % 10 == 0 | time == 0){
                output.outputPhenotypeLocation(individuals, path, time - burnIn);
                output.outputQTLGenomes(individuals, path, time - burnIn);
                output.outputNeutralGenomes(individuals, path, time - burnIn);
            }

            landscape.setIndividuals(nextIndividuals);
            offspringToAdult();

            time++;
        }

        output.outputVariances(variances, this.path);
        output.outputGenVariances(genVars, path);
    }

   //starts a simulation object
   public void runSimulationIBD() throws InterruptedException{
    //the main simulation loop will go here
    System.out.println("IBD simulation " + seed + " begins.");
    createAndShowGUI();
    addPopulation(popSize);

    int time = 0;
    int burnIn = 10000;
    double[] variances = new double[this.simLen];
    double[] genVars = new double[this.simLen];

    System.out.println("Burning in...");

    //burn-in loop
    while(time < burnIn){
        drawMutations();
        randomMating(popSize);
        offspringToAdult();

        time++;
    }

    System.out.println("Post burn-in...");

    while(time < (burnIn + simLen)){
        variances[time - burnIn] = stats.calculatePhenotypicVariance(individuals);
        genVars[time - burnIn] = stats.calculateGeneticVariance(individuals);
        
        drawMutations();
        spatialMating(this.popSize, this.individuals);

        if(time % 10 == 0 | time == 0){
            output.outputPhenotypeLocation(individuals, path, time - burnIn);
            output.outputQTLGenomes(individuals, path, time - burnIn);
            output.outputNeutralGenomes(individuals, path, time - burnIn);
        }

        landscape.setIndividuals(nextIndividuals);
        offspringToAdult();

        time++;
    }

    output.outputVariances(variances, this.path);
    output.outputGenVariances(genVars, path);
}

    //starts a simulation object
   public void runSimulationIBPXIBD() throws InterruptedException{
    //the main simulation loop will go here
    System.out.println("IBPXIBD simulation " + seed + " begins.");
    createAndShowGUI();
    addPopulation(popSize);

    int time = 0;
    int burnIn = 10000;
    double[] variances = new double[this.simLen];
    double[] genVars = new double[this.simLen];

    System.out.println("Burning in...");

    //burn-in loop
    while(time < burnIn){
        drawMutations();
        randomMating(popSize);
        offspringToAdult();

        time++;
    }

    System.out.println("Post burn-in...");

    while(time < (burnIn + simLen)){
        variances[time - burnIn] = stats.calculatePhenotypicVariance(individuals);
        genVars[time - burnIn] = stats.calculateGeneticVariance(individuals);
        
        drawMutations();
        spatialAssortativeMating();

        if(time % 10 == 0 | time == 0){
            output.outputPhenotypeLocation(individuals, path, time - burnIn);
            output.outputQTLGenomes(individuals, path, time - burnIn);
            output.outputNeutralGenomes(individuals, path, time - burnIn);
        }

        landscape.setIndividuals(nextIndividuals);
        offspringToAdult();

        time++;
    }

    output.outputVariances(variances, this.path);
    output.outputGenVariances(genVars, path);
}

    //constructor of Simulation objects that takes in all relevant parameters
    public Simulation(int scheme, int simLen, int popSize, double mutRate, double assortWidth, double disperseWidth, double spatialWidth){
        this.seed = (int) (Math.random() * 1000000);
        this.scheme = scheme;
        this.simLen = simLen;
        this.popSize = popSize;
        this.mutRate = mutRate;
        this.assortWidth = assortWidth;
        this.disperseWidth = disperseWidth;
        this.spatialWidth = spatialWidth;

        switch (scheme) {
            case 1:
                this.path = "/Users/jameson/Documents/Simulation Outputs/FullAssumptions/PAN/";
                break;
            case 2:
                this.path = "/Users/jameson/Documents/Simulation Outputs/FullAssumptions/IBP/";
                break;
            case 3:
                this.path = "/Users/jameson/Documents/Simulation Outputs/FullAssumptions/IBD/";
                break;
            case 4:
                this.path = "/Users/jameson/Documents/Simulation Outputs/FullAssumptions/IBPXIBD/";
                break;
            default:
                break;
        }
    }
}
