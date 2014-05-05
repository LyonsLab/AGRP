
import java.util.*;
import java.io.*;

public class TestCountValues{
	public static void main( String[] args) throws Exception{
        String fn = "outputFiles1/ancestorGenome_8400_9050_10997_19515.txt";
        String fnNoSubInfo = "outputFiles/ancestorGenome_8400_9050_10997_19515.txt";
        
        String[] ancestor= new String[212];
        String[] ancestorNoSubgenomeInfo = new String[171];
        
		FileReader fr1 = new FileReader(fn);
		BufferedReader br1 = new BufferedReader(fr1);
		
		String aline ="";
        for(int i = 0; i< ancestor.length; i++){
            aline = br1.readLine();
            ancestor[i] = br1.readLine();
        }
        
        fr1 = new FileReader(fnNoSubInfo);
		br1 = new BufferedReader(fr1);
		
        for(int i = 0; i< ancestorNoSubgenomeInfo.length; i++){
            aline = br1.readLine();
            ancestorNoSubgenomeInfo[i] = br1.readLine();
        }

		CommonGene cg = new CommonGene();
		GenomeInString[] twoGenomes = new GenomeInString[2];
		
        twoGenomes[0] = new GenomeInString(ancestor);
        twoGenomes[1] = new GenomeInString(ancestorNoSubgenomeInfo);
        
       
			cg.getCommonGeneForPolyploid(twoGenomes);
			BPGDistance abpg = new BPGDistance(cg.newGenomes[0].chrs, cg.newGenomes[1].chrs);
			abpg.getValue();
			System.out.println( abpg.distance + "\t"+ abpg.geneNumber);
		
	}
	
				
}
									   
	