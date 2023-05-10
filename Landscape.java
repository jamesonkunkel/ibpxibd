import java.util.ArrayList;
import javax.swing.JPanel;
import javax.swing.BorderFactory;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

public class Landscape extends JPanel {
    
    private double dimX;
    private double dimY;

    private ArrayList<Individual> individuals = new ArrayList<Individual>();
    
    //get the max X dimension of the Landscape
    public double getDimX(){
        return this.dimX;
    }

    //get the max Y dimension of the Landscape
    public double getDimY(){
        return this.dimY;
    }

    //calculate distance between any two individuals
    public double getDistance(Individual ind1, Individual ind2){
        double deltaXSquared = Math.pow((ind2.getXCoord() - ind1.getXCoord()), 2); 
        double deltaYSquared = Math.pow((ind2.getYCoord() - ind1.getYCoord()), 2);

        double distance = Math.sqrt(deltaXSquared + deltaYSquared);

        return distance;
    }

    //display locations and distance between two individuals
    public void printDistance(Individual ind1, Individual ind2){
        System.out.println("Individual 1 has X Coordinate: " + ind1.getXCoord() + " and Y Coordinate: " + ind1.getYCoord());
        System.out.println("Individual 2 has X Coordinate: " + ind2.getXCoord() + " and Y Coordinate: " + ind2.getYCoord());
        System.out.println("The distance between them is: " + getDistance(ind1, ind2));
    }

    //set the individuals present in the Landscape
    public void setIndividuals(ArrayList<Individual> individuals){
        this.individuals.clear();

        for(int i = 0; i < individuals.size(); i++){
            this.individuals.add(individuals.get(i));
        }

        repaint();
    }

    public Dimension getPreferredSize() {
        return new Dimension(500, 500);
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);       

        for(int i = 0; i < individuals.size(); i++){
            individuals.get(i).paintSquare(g);
        }
    }

    //construct new Landscape with max X and Y dimensions and an ArrayList of individuals to populate it
    public Landscape(double dimX, double dimY, ArrayList<Individual> inds){
        setBorder(BorderFactory.createLineBorder(Color.black));

        this.dimX = dimX;
        this.dimY = dimY;
        this.individuals = inds;
    }
}
