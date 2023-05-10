import java.util.Scanner;

/*
    Description: Main entry point when conducting simulations. Parameterizes simulations and then runs
    the number of replications specified. During each replication, a new Simulation object is spawned.
    The Simulation object completes it's simulation including outputs and then the next replicate begins.
*/
public class App {

    static int seed;
    static int reps;
    static int simLen;
    static int scheme;
    static int popSize;
    static double mutRate;
    static double assortWidth;
    static double disperseWidth;
    static double spatialWidth;

    //gets parameters and number of replicates for a run
    public static void input(){
        //Scanner object for System input
        Scanner scanner = new Scanner(System.in);

        //Get number of simulation replicates
        System.out.print("Replicates: ");
        reps = scanner.nextInt();

        //Get number of generations
        System.out.print("Generations: ");
        simLen = scanner.nextInt();

        System.out.print("Scheme (1=PAN, 2=IBP, 3=IBD, 4=IBPXIBD): ");
        scheme = scanner.nextInt();

        //Get population size
        System.out.print("Population size: ");
        popSize = scanner.nextInt();

        //Get mutation rate
        System.out.print("Mutation rate: ");
        mutRate = scanner.nextDouble();

        //Get width of assortative mating function
        System.out.print("Width of assortative mating function: ");
        assortWidth = scanner.nextDouble();

        //Get width of dispersal function
        System.out.print("Width of dispersal function: ");
        disperseWidth = scanner.nextDouble();

        //Get width of sptial mating function
        System.out.print("Width of spatial mating function: ");
        spatialWidth = scanner.nextDouble();

        scanner.close();
    }

    public static void main(String[] args) throws Exception {
        input();

        Simulation sim;

        //replication loop that generates new Simulation objects for every replication
        for(int i = 0; i < reps; i++){
            sim = new Simulation(scheme, simLen, popSize, mutRate, assortWidth, disperseWidth, spatialWidth);
            sim.makeOutputDirs();
            
            switch (scheme) {
                case 1:
                    sim.runSimulationPAN();
                    break;
                case 2: 
                    sim.runSimulationIBP();
                    break;
                case 3: 
                    sim.runSimulationIBD();
                    break;
                case 4: 
                    sim.runSimulationIBPXIBD();
                default:
                    break;
            }       
        }
    }
}
