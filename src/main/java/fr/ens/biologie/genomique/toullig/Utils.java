package fr.ens.biologie.genomique.toullig;

import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;

import java.io.BufferedWriter;
import java.io.IOException;

import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;

/**
 * Created by birer on 27/03/17.
 */
public class Utils {

    public Utils(){

    }
    //
    //  utils
    //



    /**
     * Method of the class TrimFastq to reverse sequence.
     * @param sequence, a sequence
     * @return a string of a reversed sequence
     */
    public static final String reverse(final String sequence) {

        if (sequence == null) {
            return null;
        }

        final char[] array = sequence.toCharArray();
        final int len = array.length;
        final StringBuilder sb = new StringBuilder(len);

        for (int i = len - 1; i >= 0; i--) {
            sb.append(array[i]);
        }
        return sb.toString();
    }

    /**
     * Get the sequence as the complement. This method work only with
     * A,T,G and C bases.
     * @param sequence sequence to reverse complement
     * @param alphabet alphabet of the sequence to reverse complement
     * @return the reverse complement sequence
     */
    public static final String complement(final String sequence,
                                           final Alphabet alphabet) {

        if (sequence == null || alphabet == null) {
            return null;
        }

        String s =
                reverseComplement(sequence,alphabet);
        final char[] array = s.toCharArray();
        final int len = array.length;
        final StringBuilder sb = new StringBuilder(len);

        for (int i = len-1; i >=  0; i--) {
            sb.append(array[i]);
        }
        return sb.toString();
    }

    //
    // Write
    //

    /**
     * Method of the class TrimFastq to write a sequence to the fasta format.
     * @param sequence, the sequence of the read
     * @param ID, the ID of the read
     * @param fastaFile, the file to write fasta sequence
     * @throws IOException
     */
    public void writeFasta(String sequence, String ID, BufferedWriter fastaFile) throws IOException {

        try {
            fastaFile.write(">"+ID+"\n");
            for (int i = 0; i <= sequence.length(); i = i + 60) {
                if (i + 60 >= sequence.length()) {
                    fastaFile.write(sequence.substring(i, sequence.length()) + "\n");
                    break;
                }
                fastaFile.write(sequence.substring(i, i + 60) + "\n");
            }
        }catch (IOException e){
            e.printStackTrace();
        }
    }

    /**
     * Method of the class TrimFastq to write a sequence to the fastq format.
     * @param ID, the ID of the read
     * @param sequence, the sequence of the read
     * @param score, the score of the read
     * @param fastqTrimBufferedWritter, the file to write fastq sequence
     * @throws IOException
     */
    public void writeFastq(String ID, String sequence, String score, BufferedWriter fastqTrimBufferedWritter)throws IOException{

        fastqTrimBufferedWritter.write("@"+ID+"\n");
        fastqTrimBufferedWritter.write(sequence+"\n");
        fastqTrimBufferedWritter.write("+\n");
        fastqTrimBufferedWritter.write(score+"\n");
    }
}
