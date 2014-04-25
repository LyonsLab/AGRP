import java.util.*;
import java.io.*;
public class EdgeForMWM{
	GenomeInString[] leaves;
	String geneContent;
    
	int[][] edgeMatrix;
	String[] geneList;
	String[] nodeString;
	
	public EdgeForMWM(GenomeInString[] l, String root){
		leaves = new GenomeInString[l.length];
		leaves = l;
		geneContent = root;
        geneList = splitBS(geneContent);
		nodeString = new String[2*geneList.length];
		for(int i = 0; i< geneList.length; i++){
			nodeString[2*i] = geneList[i]+"t";
			nodeString[2*i+1] = geneList[i]+"h";
		}
      /*  for(int i = 0; i< nodeString.length; i++){
            System.out.println(i+"\t"+nodeString[i]);
        }*/
	}
		
    public byte[][] getASetOfEdge(String[] chrs){
		byte[][] edgeMatrixHere = new byte[nodeString.length][nodeString.length];
		for(int i = 0; i<chrs.length; i++){
			String[] genes = splitBS(chrs[i]);
           // int firstNode = getFirstNode(genes[0]);
            int preNode = getSecondNode(genes[0]);
            for(int j = 1; j< genes.length; j++){
                int nodeNow = getFirstNode(genes[j]);
                edgeMatrixHere[preNode][nodeNow]++;
                edgeMatrixHere[nodeNow][preNode]++;
                preNode = getSecondNode(genes[j]);
            }
        }
       /* for(int i = 0; i< edgeMatrixHere.length; i++){
            for(int j = i+1; j<edgeMatrixHere.length; j++){
                if(edgeMatrixHere[i][j] > 0){
                    System.out.println(i+"\t"+j+"\t"+nodeString[i]+"\t"+nodeString[j]+"\t"+edgeMatrixHere[i][j]);
                }
            }
        }*/
        return edgeMatrixHere;
	}
    
    public int getFirstNode(String ag){
        String gn = ag; String  sign= "t";
        if(gn.substring(0,1).equals("-")){
            gn = gn.substring(1);
            sign = "h";
        }
        for(int i = 0; i< nodeString.length; i++){
            if(nodeString[i].equals(gn+sign)){
                return i;
            }
        }
        return -1;
    }
    
    public int getSecondNode(String ag){
        String gn = ag; String  sign= "h";
        if(gn.substring(0,1).equals("-")){
            gn = gn.substring(1);
            sign = "t";
        }
        for(int i = 0; i< nodeString.length; i++){
            if(nodeString[i].equals(gn+sign)){
                return i;
            }
        }
        return -1;
        
    }
	
	public String[] splitBS(String s){
        String[] tmp = s.trim().split(" ");
        int index = 0;
        for(int i = 0; i < tmp.length; i++){
            if(tmp[i].trim().equals("") ==false){
                index++;
            }
        }
        String[] result = new String[index];
        index = 0;
        for(int j = 0; j <tmp.length; j++){
            if(tmp[j].trim().equals("") == false){
                result[index] = tmp[j].trim();
                index++;
            }
        }
        return result;
    }
	
	
}