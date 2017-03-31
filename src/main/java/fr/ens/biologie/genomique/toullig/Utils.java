package fr.ens.biologie.genomique.toullig;

import fr.ens.biologie.genomique.eoulsan.bio.Alphabet;

import static fr.ens.biologie.genomique.eoulsan.bio.Sequence.reverseComplement;

/**
 * Class of Utils methods. Created by birer on 27/03/17.
 */
public class Utils {

  public Utils() {

  }
  //
  // utils
  //

  /**
   * Method of the class Utils to reverse sequence.
   * @param sequence, a sequence
   * @return a string of a reversed sequence
   */
  public static String reverse(final String sequence) {

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
   * Get the sequence as the complement. This method work only with A,T,G and C
   * bases.
   * @param sequence sequence to reverse complement
   * @param alphabet alphabet of the sequence to reverse complement
   * @return the reverse complement sequence
   */
  public static String complement(final String sequence,
      final Alphabet alphabet) {

    if (sequence == null || alphabet == null) {
      return null;
    }

    String s = reverseComplement(sequence, alphabet);
    final char[] array = s.toCharArray();
    final int len = array.length;
    final StringBuilder sb = new StringBuilder(len);

    for (int i = len - 1; i >= 0; i--) {
      sb.append(array[i]);
    }
    return sb.toString();
  }
}
