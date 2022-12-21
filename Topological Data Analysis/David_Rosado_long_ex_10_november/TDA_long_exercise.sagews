︠3c7d1076-1506-4ed4-9494-bf79be2b689es︠
#Let us define the simplicial complex K and L

K = SimplicialComplex([[1,2,4],[1,2,5],[1,3,5],[1,3,6],[1,4,6],[2,3,4],[2,3,6],[2,5,6],[3,4,5],[4,5,6]]) #define the simplicial complex
L = SimplicialComplex([[0,1,4],[0,1,5],[0,2,3],[0,2,7],[0,3,5],[0,4,7],[1,2,6],[1,2,8],[1,4,8],[1,5,6],[2,3,6],[2,7,8],[3,4,6],[3,4,8],[3,5,8],[4,6,7],[5,6,7],[5,7,8]])

#Let us check the maximal faces

K.facets()#maximal faces of the simplicial complex
L.facets()#maximal faces of the simplicial complex

︡6c98eecc-b203-4077-a884-8ab15c2f8ef4︡{"stdout":"{(1, 3, 5), (2, 3, 4), (1, 2, 4), (1, 3, 6), (2, 3, 6), (1, 2, 5), (2, 5, 6), (3, 4, 5), (4, 5, 6), (1, 4, 6)}\n"}︡{"stdout":"{(2, 7, 8), (2, 3, 6), (0, 3, 5), (1, 4, 8), (0, 4, 7), (0, 1, 4), (4, 6, 7), (0, 2, 3), (0, 2, 7), (3, 5, 8), (5, 7, 8), (3, 4, 6), (0, 1, 5), (5, 6, 7), (1, 2, 8), (3, 4, 8), (1, 2, 6), (1, 5, 6)}\n"}︡{"done":true}
︠82f35cd4-b0d2-4f85-885e-6709837d0879s︠
K.homology(None,ZZ,None,False,False,'pari',False,False)#homology group with coefficients in Z
print('\n')
L.homology(None,ZZ,None,False,False,'pari',False,False)#homology group with coefficients in Z
︡47720668-21af-4e8e-9f64-4f4882f85390︡{"stdout":"{0: Z, 1: C2, 2: 0}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Z, 1: Z x Z, 2: Z}\n"}︡{"done":true}
︠28e3d354-7a91-4e3d-bc69-e66ff92fb2dfs︠
K.homology(None,QQ,None,False,False,'pari',False,False)#homology group with coefficients in Q
print('\n')
L.homology(None,QQ,None,False,False,'pari',False,False)#homology group with coefficients in Q
︡2cb2e1e7-4db8-4905-81aa-6560e647e7b7︡{"stdout":"{0: Vector space of dimension 1 over Rational Field, 1: Vector space of dimension 0 over Rational Field, 2: Vector space of dimension 0 over Rational Field}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Vector space of dimension 1 over Rational Field, 1: Vector space of dimension 2 over Rational Field, 2: Vector space of dimension 1 over Rational Field}\n"}︡{"done":true}
︠12116d4c-8e01-4a3b-a16e-19a709b65445s︠
K.homology(None,GF(2),None,False,False,'pari',False,False)#homology group with coefficients in GF(2)
print('\n')
L.homology(None,GF(2),None,False,False,'pari',False,False)#homology group with coefficients in GF(2)
︡19229003-3b8a-4e5d-8b78-9e4d5c8fcc89︡{"stdout":"{0: Vector space of dimension 1 over Finite Field of size 2, 1: Vector space of dimension 1 over Finite Field of size 2, 2: Vector space of dimension 1 over Finite Field of size 2}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Vector space of dimension 1 over Finite Field of size 2, 1: Vector space of dimension 2 over Finite Field of size 2, 2: Vector space of dimension 1 over Finite Field of size 2}\n"}︡{"done":true}
︠71e749e6-e7af-4431-9603-02d4cfdd8c5es︠
#Let us compute the homology of the real projective plane and the torus
RP = simplicial_complexes.RealProjectivePlane()
T = simplicial_complexes.Torus()
︡ed5fae43-56ca-4a33-adb3-0ad6a5ec583a︡{"done":true}
︠1ef431bc-d16b-4f1a-bea3-2582ff0aa954s︠
RP.homology(None,ZZ,None,False,False,'pari',False,False)#homology group with coefficients in Z
print('\n')
T.homology(None,ZZ,None,False,False,'pari',False,False)#homology group with coefficients in Z
︡3de8f817-a5c8-4446-af39-46ffb8582ad5︡{"stdout":"{0: Z, 1: C2, 2: 0}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Z, 1: Z x Z, 2: Z}\n"}︡{"done":true}
︠96989156-f6c7-494b-ac90-71c040ab394bs︠
RP.homology(None,QQ,None,False,False,'pari',False,False)#homology group with coefficients in Q
print('\n')
T.homology(None,QQ,None,False,False,'pari',False,False)#homology group with coefficients in Q
︡466aa538-8814-4806-a82b-364d40283185︡{"stdout":"{0: Vector space of dimension 1 over Rational Field, 1: Vector space of dimension 0 over Rational Field, 2: Vector space of dimension 0 over Rational Field}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Vector space of dimension 1 over Rational Field, 1: Vector space of dimension 2 over Rational Field, 2: Vector space of dimension 1 over Rational Field}\n"}︡{"done":true}
︠9ea39ff7-55ff-41ad-bb89-fd12fe416991s︠
RP.homology(None,GF(2),None,False,False,'pari',False,False)#homology group with coefficients in GF(2)
print('\n')
T.homology(None,GF(2),None,False,False,'pari',False,False)#homology group with coefficients in GF(2)
︡692903bb-e61a-48d8-be81-09b6c95971f2︡{"stdout":"{0: Vector space of dimension 1 over Finite Field of size 2, 1: Vector space of dimension 1 over Finite Field of size 2, 2: Vector space of dimension 1 over Finite Field of size 2}\n"}︡{"stdout":"\n\n"}︡{"stdout":"{0: Vector space of dimension 1 over Finite Field of size 2, 1: Vector space of dimension 2 over Finite Field of size 2, 2: Vector space of dimension 1 over Finite Field of size 2}\n"}︡{"done":true}
︠d7b7d0d2-cba5-4bf9-b5f7-749a9597925f︠









