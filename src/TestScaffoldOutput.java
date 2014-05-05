
import java.util.*;
import java.io.*;

public class TestScaffoldOutput{
	public static void main( String[] args) throws Exception{
		
        if (args.length < 8) {
            System.out.println(Usage());
            System.exit(1);
        }
        
        String[] genomeIndex = null;
        String[] weights = null;
        String[] ploidyArray = null;
        int weightToInclude = 1;
        String minimumLengthString = null;
        String inputGenomes = null;
        String output_directory = ".";
        
        // Input files
        File [] scaffoldFiles = null;
        File [] binFiles = null;
        String contigs2genes = null;
        String outputFileName = "ancestorGenome.txt";
        
        for(int i = 0; i < args.length; i += 2) {
            String argument = args[i].substring(1);
            String option = args[i + 1];
            
            switch(argument) {
                case "im":
                    File dir1 = new File(option);
                    scaffoldFiles = dir1.listFiles();
                    for(int j = 0; j < scaffoldFiles; j++) {
                        System.out.println(scaffoldFiles[j].getName());
                        
                    }
                    break;
                case "ib":
                    File dir2 = new File(option);
                    binFiles = dir2.listFiles();
                    
                    for(int j = 0; j < binFiles; j++) {
                        System.out.println(binFiles[i].getName());
                    }
                    break;
                case "cg":
                    contigs2genes = option;
                    break;
                case "o":
                    outputFileName = option;
                    break;
                default:
                    System.out.println("Unknown argument");
                    System.out.println(Usage());
                    System.exit(1);
            }
        }
        
        if (contigs2genes == null || scaffoldFiles == null || binFiles == null) {
            System.out.println("Please specify all the command line arguments");
            System.out.println(Usage());
            System.exit(1);
        }
        
        int[] allGenomeIndex = new int[4];
        allGenomeIndex[0] = 8400; allGenomeIndex[1] = 9050; allGenomeIndex[2] = 10997; allGenomeIndex[3] = 19515;
        // read in contig
        FileReader fr = new FileReader(contigs2genes);
        BufferedReader br = new BufferedReader(fr);
        String aline = br.readLine();
        String[] info = aline.split("\t");
        int totalContig = Integer.parseInt(info[1]);
        String[] contigs = new String[totalContig];
        for(int g = 0; g< contigs.length; g++){
            aline = br.readLine();
            contigs[g] = br.readLine();
        }
        
       int ancChrNumber = binFiles.length;
        
        GenomeInString[] ancestorChrs = new GenomeInString[ancChrNumber];
        int ancTotalFragment = 0;
        for(int c = 1; c< ancChrNumber+1; c++){
            // read in gene content
          //  String geneContentFile = "outputFiles/leaveGenomesInContigForAAncChr"+new Integer(c).toString()+".txt";
            String geneContent = "";
           // System.out.println("geneContent file "+ geneContentFile);
            fr = new FileReader(binFiles[c]);
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
         //   String inputFile = "outputFiles/scaffoldOutput"+new Integer(c).toString()+".txt";
            fr = new FileReader(scaffoldFiles[c]);
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
            String[] aChr = scaffoldFinal.getGenomeHighEdgeWeight(0);
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
      //  String finalOutputFile = "outputFiles/ancestorGenome";
       /* for(int i =0; i<  allGenomeIndex.length; i++){
            finalOutputFile = finalOutputFile+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        finalOutputFile = finalOutputFile+".txt";*/
        FileWriter fstream = new FileWriter(outputFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        for(int i = 0; i< finalAnc.length; i++){
            fbw.write("chr "+ i+"\n"+finalAnc[i]+"\n");
		}
        
        fbw.flush();
        fbw.close();
        
	}
    
    public static String Usage() {
        return "TODO";
    }
}