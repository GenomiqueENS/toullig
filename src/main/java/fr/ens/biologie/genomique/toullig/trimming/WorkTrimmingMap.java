// package fr.ens.biologie.genomique.toullig.trimming;
//
//
// import java.util.Objects;
//
/// **
// * Class of the Working trimming Map
// * Created by birer on 03/04/17.
// */
// public class WorkTrimmingMap {
//
//
// static class InfomationRead {
// private String sequence;
// private String quality;
//
//
//
//
// InfomationRead(Integer id, String name ) {
// this.id = id;
// this.name = name;
//
//
//
//
//
//
// }
//
// @Override
// public int hashCode() {
// // this ensures all hashcodes are between 0 and 15
// return Objects.hash(id,name );
// }
//
// @Override
// public boolean equals(Object obj) {
//
// if(obj==null){
// return false;
// }
// if(!(obj instanceof InfomationRead)){
// return false;
// }
//
// InfomationRead otherEmp = (InfomationRead) obj;
// return this.name.equals(otherEmp.name);
// }
//
// @Override
// public String toString() {
// return this.id + "-" + name;
// }
// }
// }
