/*
 * Parallel merge sort algorithm implemented in Java
 */
package lab1;

import java.util.Scanner;
import java.util.InputMismatchException;
import java.util.ArrayList;
import java.util.Random;

class run {

    public static void main(String[] args) {
        tester t = new tester();
        t.start();    
    }  
}

/*
* A Tester class to test this merge sort program
*/
class tester {
    //public int myArr[] = {5, 160, 2, 87, 299, 150, 999};
    public int size;
    public ArrayList<Integer> myArr = new ArrayList<>();
    tester(){
        System.out.println("LAB1 - array merge sort ...");
        // initialise the array list
        this.arr_init();
    }
    
    public void arr_init() {
        int nThread = 0;
        long randSeed = 0;
        int arrSize = 0;
        
        System.out.println("Pleae enter the size of array, a rand seed: ");
        // get array size from user inputs
        Scanner in = new Scanner(System.in);
        try {
                arrSize = in.nextInt();
                randSeed = in.nextLong();
                in.close();
        }
        catch (InputMismatchException e) {
            throw new InputMismatchException("Please enter an ");
        }
        
        Random rand = new Random(randSeed);
        
        // initialise the array
        for (int i=0; i<arrSize; i++) {
            this.myArr.add(rand.nextInt(1000 - 1) + 1);
        }
        
        // print out original array
        for(int i = 0; i < arrSize; i++) {  
            if (i == (arrSize - 1)) {
                System.out.print(this.myArr.get(i) + " ");
            }
            else {
                System.out.print(this.myArr.get(i) + ", ");
            }
        }
        System.out.println();
    }
    
    public void start() {
        // instentiate an object of the sorter class and apply the sort
        int length = myArr.size();
        int rIndex = length - 1;
        sorter st = new sorter(this.myArr, 0, rIndex);
        
        System.out.println("\n === After sorting ===\n");
        
        // print the sorted array
        for(int i = 0; i < this.myArr.size(); i++) {  
            if (i == (this.myArr.size() - 1)) {
                System.out.print(this.myArr.get(i) + " ");
            }
            else {
                System.out.print(this.myArr.get(i) + ", ");
            }
        }
        System.out.println();
    }    
}




 

