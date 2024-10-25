with Symmetry_Operations; use Symmetry_Operations;

package Shmueli_Matrices is
   
   -- This package describes the rotation matrices necessary for
   --  interpretation of Shmueli explicit space group symbols [1,2].
   
   -- Refs.:
   
   -- [1] Shmueli, U. (1984) Space-group algorithms. I. The space
   --  group and its symmetry elements. Acta Crystallographica Section
   --  A Foundations of Crystallography 40(5), 559-567. International
   --  Union of Crystallography (IUCr). DOI:
   --  https://doi.org/10.1107/s0108767384001161
   --
   -- Table 2. Point-group generators (proper rotations)
   
   -- [2] U. Shmueli (ed.) (2001) International Tables for
   --  Crystallography, Vol. B. Reciprocal Space. (c) International
   --  Union of Crystallography 2001. Published by Kluwer Academic
   --  Publishers, Dordrecht/Boston/London. p. 108.
   
   subtype Crystallographic_Integer is Integer range -1 .. 1;
   
   type Rotation_Matrix is array (1 .. 3, 1 .. 3) of Crystallographic_Integer;
   
   M_1A : Rotation_Matrix :=
     (
      (1 => 1, others => 0),
      (2 => 1, others => 0),
      (3 => 1, others => 0)
     );
   
   M_2A : Rotation_Matrix :=
     (
      (1 =>  1, others => 0),
      (2 => -1, others => 0),
      (3 => -1, others => 0)
     );
   
   M_2B : Rotation_Matrix :=
     (
      (1 => -1, others => 0),
      (2 =>  1, others => 0),
      (3 => -1, others => 0)
     );
   
   M_2C : Rotation_Matrix :=
     (
      (1 => -1, others => 0),
      (2 => -1, others => 0),
      (3 =>  1, others => 0)
     );
   
   M_2D : Rotation_Matrix :=
     (
      ( 0,  1,  0),
      ( 1,  0,  0),
      ( 0,  0, -1)
     );
   
   M_2E : Rotation_Matrix :=
     (
      ( 0, -1,  0),
      (-1,  0,  0),
      ( 0,  0, -1)
     );
   
   -- "Note that the axes 2F, 2G, 3C and 6C appear in trigonal and
   --  hexagonal space groups":
   
   M_2F : Rotation_Matrix :=
     (
      ( 1, -1,  0),
      ( 0, -1,  0),
      ( 0,  0, -1)
     );
   
   M_2G : Rotation_Matrix :=
     (
      ( 1,  0,  0),
      ( 1, -1,  0),
      ( 0,  0, -1)
     );
   
   M_3Q : Rotation_Matrix :=
     (
      ( 0,  0,  1),
      ( 1,  0,  0),
      ( 0,  1,  0)
     );
   
   M_3C : Rotation_Matrix :=
     (
      ( 0, -1,  0),
      ( 1, -1,  0),
      ( 0,  0,  1)
     );
   
   M_4C : Rotation_Matrix :=
     (
      ( 0, -1,  0),
      ( 1,  0,  0),
      ( 0,  0,  1)
     );
   
   M_6C : Rotation_Matrix :=
     (
      ( 1, -1,  0),
      ( 1,  0,  0),
      ( 0,  0,  1)
     );
   
   procedure Decode_Shmueli_Symbol_Translation
     (
      Translation_Code : in String;
      Translation      : out Crystallographic_Translation
     );
   
end Shmueli_Matrices;
