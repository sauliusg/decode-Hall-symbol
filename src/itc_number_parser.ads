with Symmetry_Operations; use Symmetry_Operations;

package ITC_Number_Parser is
   
   ITC_NUMBER_SYMBOL_NOT_FOUND : exception;

   Debug_Print_Matrices : Boolean := False;
   
   subtype Space_Group_Number_Type is Positive range 1 .. 230;
   
   function Lookup_ITC_Number (SG_Number : in Space_Group_Number_Type)
                              return Symmetry_Operator_Array;
   
   function Parse_ITC_Number (ITC_Number_And_Cell : in String)
                              return Symmetry_Operator_Array;
   
end ITC_Number_Parser;
