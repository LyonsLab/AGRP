
import java.util.*;
import java.io.*;

public class TestScaffoldOutputNoSubgenomeInfo{
	public static void main( String[] args) throws Exception{
		
        int[] allGenomeIndex = new int[4];
        allGenomeIndex[0] = 8400; allGenomeIndex[1] = 9050; allGenomeIndex[2] = 10997; allGenomeIndex[3] = 19515;
        // read in contig
        String contigFile = "outputFiles/contig_8400_9050_10997_19515.txt";
        FileReader fr = new FileReader(contigFile);
        BufferedReader br = new BufferedReader(fr);
        String aline = br.readLine();
        String[] info = aline.split("\t");
        int totalContig = Integer.parseInt(info[1]);
        String[] contigs = new String[totalContig];
        for(int g = 0; g< contigs.length; g++){
            aline = br.readLine();
            contigs[g] = br.readLine();
        }
        
        int ancChrNumber = 1;
        
        GenomeInString[] ancestorChrs = new GenomeInString[ancChrNumber];
        int ancTotalFragment = 0;
        for(int c = 1; c< ancChrNumber+1; c++){
            // read in gene content
            String geneContentFile = "outputFiles/leaveGenomesInContig.txt";
            String geneContent = "";
            System.out.println("geneContent file "+ geneContentFile);
            fr = new FileReader(geneContentFile);
            br = new BufferedReader(fr);
            aline = br.readLine();
            aline = br.readLine();
            geneContent = br.readLine();
            System.out.println(" gene content " + geneContent);
            int geneContentLength = geneContent.length();
            System.out.println("geneContenet length\t"+geneContentLength);
            // read in adj
            String[] edges = new String[geneContent.length()/3];
            int index = 0;
            String inputFile = "outputFiles/scaffoldOutput.txt";
            fr = new FileReader(inputFile);
            br = new BufferedReader(fr);
            aline = br.readLine();
            while(aline!=null && aline.substring(0,5).equals("total")==false){
                edges[index] = aline;
                index++;
                aline = br.readLine();
            }
            System.out.println("total edges "+ index);
            System.out.println("---------"+aline);
            

            GenomeAdj scaffoldFinal = new GenomeAdj(geneContent);
            scaffoldFinal.initialValue(edges);
            int weightToInclude = 2;
            String[] aChr = scaffoldFinal.getGenomeHighEdgeWeight(weightToInclude);
            ancTotalFragment= ancTotalFragment+ aChr.length;
            String[] chrInGene = scaffoldFinal.getAncestorInGene(aChr, contigs);
           // for(int i = 0; i< chrInGene.length; i++){
           //     System.out.println("chr "+i+"\n"+chrInGene[i]);
           // }
            ancestorChrs[c-1] = new GenomeInString(chrInGene);
        }
        //
        String[] finalAnc = new String[ancTotalFragment];
        int indexhere = 0;
        for(int i = 0; i< ancestorChrs.length;i++){
            for(int j = 0; j< ancestorChrs[i].chrs.length;j++){
                finalAnc[indexhere] = ancestorChrs[i].chrs[j];
                indexhere++;
            }
            
        }
        
        GenomeInString ancestor = new GenomeInString(finalAnc);
       int totalGeneInAnc =  ancestor.countGeneNumber();
        System.out.println("total gene in ancestor\t"+totalGeneInAnc);
        
        // final output
        String finalOutputFile = "outputFiles/ancestorGenome";
        for(int i =0; i<  allGenomeIndex.length; i++){
            finalOutputFile = finalOutputFile+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        finalOutputFile = finalOutputFile+".txt";
        FileWriter fstream = new FileWriter(finalOutputFile,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        for(int i = 0; i< finalAnc.length; i++){
            fbw.write("chr "+ i+"\n"+finalAnc[i]+"\n");
		}
        
        fbw.flush();
        fbw.close();
        
		
        
	}
}