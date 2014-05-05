
import java.util.*;
import java.io.*;

public class GeneInfo{
	String name;
	String newName;
	int chrNumber;
	int positionStart;
	int positionEnd;
	String orientation;
	int genomeIndex;
	int ploidyNumber = 1;
    
	public GeneInfo(String gn, int gi){
		name = gn;
		genomeIndex = gi;
	}
    
	public GeneInfo(){
	}
	public GeneInfo(GeneInfo g){
		name = g.name;
		newName = g.newName;
		chrNumber = g.chrNumber;
		positionStart = g.positionStart;
		positionEnd = g.positionEnd;
		orientation = g.orientation;
		genomeIndex=g.genomeIndex;
		ploidyNumber = g.ploidyNumber;
	}
	
	public boolean sameGene(GeneInfo ag){
		if(genomeIndex==ag.genomeIndex){
			//if(ag.chrNumber==chrNumber){
				if(name.equals(ag.name)){
				//if(chrNumber!=ag.chrNumber || positionStart!=ag.positionStart || positionEnd!=ag.positionEnd){
				//	System.out.println("same name, different position");this.printGene();ag.printGene();
				//	System.out.print("&");
				//}
				return true;
			//}
			}}
		return false;
	}
	
	public void print(){
		System.out.println(name + "\t" + newName+  "\t" + chrNumber + "\t" + positionStart + "\t" + positionEnd+"\t" + orientation+"\t"+genomeIndex+"\t"+ploidyNumber);
	}

}
