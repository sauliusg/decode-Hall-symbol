with Symmetry_Operations; use Symmetry_Operations;

package Change_Of_Basis is
   
   -- Parse the Change of Basis (CoB) operator strings in the form
   --  "(a-b,a+b,c)" or "(x-y,x+y,z)" and return a matrix that encodes
   --  this operator. The CoB is described in [1,2].
   --
   -- The 'abc' for is a transpose inverse of the 'xyz' form.
   --
   -- Refs.:
   --
   -- [1] Zwart, P. H.; Grosse-Kunstleve, R. W.; Lebedev, A. A.;
   --  Murshudov, G. N. & Adams, P. D. (2007) Surprises and pitfalls
   --  arising from (pseudo)symmetry. Acta Crystallographica Section D
   --  Biological Crystallography 64(1), 99-107. International Union
   --  of Crystallography (IUCr). DOI:
   --  https://doi.org/10.1107/s090744490705531x
   --
   -- [2] Sydney R. Hall, Ralf W. Grosse-Kunstleve (1996) "Concise
   --  Space-Group Symbols". URL:
   --  https://cci.lbl.gov/sginfo/hall_symbols.html [accessed:
   --  2022-06-14T15:24+03:00]
   --   
   -- [3] International Tables Volume B (2010), "Symmetry in
   --  reciprocal space". Section 1.4., Appendix A1.4.2. Space-group
   --  symbols for numeric and symbolic computations, URL:
   --  https://onlinelibrary.wiley.com/iucr/itc/B/ [accessed:
   --  2022-06-14T15:35+03:00]
   
   WRONG_CHANGE_OF_BASIS : exception;
   
   procedure Get_Change_Of_Basis (
                                  Symbol : in String;
                                  Pos : in out Integer;
                                  Change_Of_Basis : out Symmetry_Operator
                                 );   
end Change_Of_Basis;
