
import java.util.*;
import java.io.*;
import java.util.regex.Pattern;

public class TestGetContigInput{
	public static void main( String[] args) throws Exception{
        if (args.length < 8) {
            System.out.println(Usage());
            System.exit(1);
        }
        
        String[] genomeIndex = null;
        String[] weightSchemeString = null;
        String inputGenomes = null;
        String weightOfAdjacency = null;
        String output_directory = ".";
        
        for(int i = 0; i < args.length; i += 2) {
            String argument = args[i].substring(1);
            String option = args[i + 1];
            
            switch(argument) {
                case "g":
                    String genomeIndexString = option;
                    genomeIndex = option.split(Pattern.quote(","));
                    break;
                case "w":
                    String weightString = option;
                    weightSchemeString = weightString.split(Pattern.quote(","));
                    break;
                case "wa":
                    weightOfAdjacency = option;
                    break;
                case "i":
                    inputGenomes = option;
                    break;
                case "o":
                    output_directory = option;
                    break;
                default:
                    System.out.println("Unknown argument");
                    System.out.println(Usage());
                    System.exit(1);
            }
        }

        if (weightSchemeString == null || genomeIndex == null || inputGenomes == null ||
            weightOfAdjacency == null) {
            System.out.println("Please specify all the command line arguments");
            System.out.println(Usage());
            System.exit(1);
        }
        
        int numberOfGenomes =genomeIndex.length;
        int weightToInclude = Integer.parseInt(weightOfAdjacency);
       // int[] allGenomeIndex = new int[numberOfGenomes];
        int[] allGenomeIndex = new int[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            allGenomeIndex[i] = Integer.parseInt(genomeIndex[i]);
        }
        int[] weightScheme = new int[weightSchemeString.length];
        for(int i = 0; i< weightSchemeString.length; i++){
            weightScheme[i] = Integer.parseInt(weightSchemeString[i]);
        }
        
        GenomeInString[] leaveGenomes = new GenomeInString[numberOfGenomes];
        
        String leaveGenomesFile = inputGenomes;
        FileReader fr1 = new FileReader(leaveGenomesFile);
        BufferedReader br1 = new BufferedReader(fr1);
        String aline1 = br1.readLine();
        String[] info = aline1.split("\t");
        int numberOfGenes =Integer.parseInt(info[1]);
        System.out.println("numberOfGenes\t"+numberOfGenes);
        int[] chrNumber = new int[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            aline1 = br1.readLine();
            System.out.println(aline1);
            info = aline1.split("\t");
            chrNumber[i]  =Integer.parseInt(info[3]);
        }
        
        for(int i = 0; i< leaveGenomes.length; i++){
            aline1 = br1.readLine();  // genome name
            String[] ag = new String[chrNumber[i]];
            for(int j = 0; j< ag.length; j++){
                aline1 = br1.readLine();
                ag[j] = br1.readLine();
            }
            leaveGenomes[i] = new GenomeInString(ag);
           
        }
        
		// get Gene Content
		String geneContent = "";
		for(int i = 1; i<numberOfGenes+1; i++){
			geneContent = geneContent+"  "+ new Integer(i).toString();
		}
	
		EdgeForMWM allEdges = new EdgeForMWM(leaveGenomes, geneContent);
        
        int[][] edgeMatrix = new int[allEdges.nodeString.length][allEdges.nodeString.length];
        
        for(int i = 0; i< leaveGenomes.length; i++){
            byte[][] edgeMatrix1 = allEdges.getASetOfEdge(leaveGenomes[i].chrs);
            for(int j = 0; j< edgeMatrix.length; j++){
                for(int k = 0; k< edgeMatrix.length; k++){
                    int weighthere = weightScheme[i]*edgeMatrix1[j][k];
                    edgeMatrix[j][k] = edgeMatrix[j][k] +weighthere;
                }
            }
        }
        
       /* for(int i = 0; i< edgeMatrix.length-1; i++){
            for(int j = i+1; j< edgeMatrix.length; j++){
                if(edgeMatrix[i][j]>0){
                    System.out.println(i+"\t"+j+"\t"+allEdges.nodeString[i]+"\t"+allEdges.nodeString[j]+"\t"+edgeMatrix[i][j]);
                }
            }
        }
     */
        
        /*int[] edgeWeight = new int[2000];
		for(int i = 0; i< edgeMatrix.length-1; i++){
			for(int j = i+1; j< edgeMatrix.length; j++){
				edgeWeight[edgeMatrix[i][j]]++;
			}
		}
	
		System.out.println("Weight\tnumberOFAdj");
		for(int i = 0; i< edgeWeight.length; i++){
			if(edgeWeight[i]!=0){
				System.out.println(i+ "\t" + edgeWeight[i]);
			}
		}
		
		int totalNumberOfEdge = 0;
		for(int i = 0; i< edgeMatrix.length-1; i++){
			for(int j = i+1; j< edgeMatrix.length; j++){
				if(edgeMatrix[i][j]!=0){
					totalNumberOfEdge++;
				}
			}
		}*/
        
        
        String outPutFileName = output_directory + "/contigInput";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
      /*  outPutFileName = outPutFileName+".py";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        String fileName = "mwmatching.py";
        fbw.write("import time\n");
        fbw.write("import sys\n");
        fbw.write("start = time.clock()\n");
        fbw.write("sys.setrecursionlimit(4000)\n");
        fr = new FileReader(fileName);
        br = new BufferedReader(fr);
        aline = br.readLine();
        while(aline!=null){
            fbw.write(aline+"\n");
            aline = br.readLine();
        }
        fbw.write("maxWeightMatching([ ");
        for(int i = 0; i< edgeMatrix.length-1; i++){
            for(int j = i+1; j< edgeMatrix.length; j++){
                if(edgeMatrix[i][j]>=weightToInclude){
                    fbw.write("("+i+","+j+","+edgeMatrix[i][j]+"),");
                }
            }
        }
        fbw.write("])\n");
        fbw.write("end = time.clock()\n");
        fbw.write("print end - start\n");
       */
        outPutFileName = outPutFileName+".txt";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        for(int i = 0; i< edgeMatrix.length-1; i++){
            for(int j = i+1; j< edgeMatrix.length; j++){
                if(edgeMatrix[i][j]>=weightToInclude){
                    fbw.write(i+"\t"+j+"\t"+edgeMatrix[i][j]+"\n");
                }
            }
        }
        fbw.flush();
        fbw.close();
    
    
	}
    
    public static String Usage() {
        return "TODO";
    }
}