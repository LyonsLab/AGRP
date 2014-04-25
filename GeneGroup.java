import java.util.*;
import java.io.*;
public class GeneGroup{
	
	int genomeIndex;
	int geneNumber;
	String[] geneName;
	int chr;
	double averagePos;
	int[] pos;
	int colorCode;
	int blockIndex;
	String sign;
	
	public GeneGroup(){}
    
    public GeneGroup(GeneGroup ag){
        genomeIndex = ag.genomeIndex;
        geneNumber = ag.geneNumber;
        geneName = new String[ag.geneName.length];
        geneName = ag.geneName;
        chr = ag.chr;
        averagePos = ag.averagePos;
        pos = new int[ag.pos.length];
        pos = ag.pos;
        colorCode = ag.colorCode;
        blockIndex = ag.blockIndex;
        sign = ag.sign;
    }
	public GeneGroup(int gi, int gn, String[] gns, int c, int[] p, int cc, int bi, String s ){
		genomeIndex = gi;
		geneNumber = gn;
		geneName = new String[gns.length];
		geneName = gns;
		chr = c;
		pos = new int[p.length];
		pos = p;
		colorCode = cc;
		blockIndex = bi;
		sign = s;
		double sum = 0;
		for(int i = 0; i< pos.length; i++){
			sum = sum+pos[i];
		}
		averagePos = sum/pos.length;
	}
	
	public void printAGeneGroup(){
		System.out.println("==="+genomeIndex+"\t"+geneNumber+"\t"+chr+"\t"+colorCode+"\t"+blockIndex+"\t\t"+ averagePos);
		for(int i = 0; i< geneName.length; i++){
			System.out.print(geneName[i]+"\t");
		}
		System.out.println();
		for(int i = 0; i< geneName.length; i++){
			System.out.print(pos[i]+"\t");
		}
		System.out.println();
		
	}
	
	
}