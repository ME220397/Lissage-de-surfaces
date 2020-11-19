#include "mainwindow.h"
#include "ui_mainwindow.h"

/* **** début de la partie à compléter **** */

double cot(double teta){
    return std::cos(teta)/std::sin(teta);
}

void scale(MyMesh * _mesh){
    Vec3f new_coords;
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        new_coords = _mesh->point(*v);
        _mesh->set_point(*v, new_coords*2);
    }
}

void unscale(MyMesh * _mesh){
    Vec3f new_coords;
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        new_coords = _mesh->point(*v);
        _mesh->set_point(*v, new_coords*0.5);
    }
}


void MainWindow::flou_de_diffusion(MyMesh *_mesh, VertexHandle v, MyMesh::Point fv, double f){
    // On calcul le nouveau point v
    MyMesh::Point new_v = _mesh->point(v) + f * fv;
    qDebug() << "old v :"  << _mesh->point(v)[0] << _mesh->point(v)[1] << _mesh->point(v)[2];
    qDebug() << " vecteur directeur : " << fv[0] << fv[1] << fv[2];
    qDebug() << "new v :"  << new_v[0] << new_v[1] << new_v[2];
    _mesh->set_point(v, new_v);
}

void MainWindow::operateur_laplace_beltrami(MyMesh* _mesh, int choix, double h, double _y){
    // on va maintenant appliqué le flou de diffusion en chaque point que l'on recupère dans un vecteur
    QVector<MyMesh::Point> fvi;
    if(choix == COTANGENTE){
        for (MyMesh::VertexIter current = _mesh->vertices_begin(); current != _mesh->vertices_end(); current++) {
            fvi.append(laplace_beltrami_cot(_mesh, *current));
            //MyMesh::Point p = laplace_beltrami_cot(_mesh, *current);
            //flou_de_diffusion(_mesh, *current, p, _y);
            //break;
        }
    }
    else{
        for (MyMesh::VertexIter current = _mesh->vertices_begin(); current != _mesh->vertices_end(); current++) {
            fvi.append(laplace_beltrami_uni(_mesh, *current));
            //MyMesh::Point p = laplace_beltrami_cot(_mesh, *current);
            //flou_de_diffusion(_mesh, *current, p, _y);
            //break;
        }
    }

    for (int i = 0; i<fvi.length(); i++) {
        // On applique le flou de diffusion à tous les points
        VertexHandle v = _mesh->vertex_handle(i);
        flou_de_diffusion(_mesh, v, fvi.at(i), _y);
    }
}

MyMesh::Point MainWindow::laplace_beltrami_uni(MyMesh* _mesh, VertexHandle v){
    // A completer
    MyMesh::Point direction;
    direction *= 0; // Initialisation
    MyMesh::Point u;
    u *= 0; // Initialisation
    MyMesh::Point somme;
    somme *= 0; // Initialisation
    int n_voisin = 0;
    for(MyMesh::VertexVertexIter vi = _mesh->vv_iter(v); vi.is_valid(); vi++){
        // 3. Calcul des du vecteur directeur vvi
        direction = calc_vector_v_vi(_mesh, v, *vi);
        u = direction;
        //u.normalize();
        somme += u;
        n_voisin++;
    }
    MyMesh::Point laplacien_cot_v = 1/(n_voisin) * somme;
    //laplacien_cot_v.normalize();

    return laplacien_cot_v;

}

MyMesh::Point MainWindow::laplace_beltrami_cot(MyMesh *_mesh, VertexHandle v){

    // L'approcimation cotengentielle de laplace beltrami peut se diviser en 3 parties.

    // 1. Calcul de l'aire
    MyMesh::Scalar aire = neighboring_faces_area(_mesh, v);

    // Pour tout les voisins de v, on va recuperer la somme des poids multiplié
    // par le vecteur directeur vvi que l'on va noter u;
    MyMesh::Scalar poids;
    MyMesh::Point direction;
    direction *= 0; // Initialisation
    MyMesh::Point u;
    u *= 0; // Initialisation
    MyMesh::Point somme;
    somme *= 0; // Initialisation
    for(MyMesh::VertexVertexIter vi = _mesh->vv_iter(v); vi.is_valid(); vi++){
        // 2. Calcul des points cot
        poids = calcul_poids_cot(_mesh, v, *vi);
        // 3. Calcul des du vecteur directeur vvi
        direction = calc_vector_v_vi(_mesh, v, *vi);
        u = poids*direction;
        //u.normalize();
        somme += u;
    }
    //somme.normalize();
    //Maintenant que nous avons tout calculé, on fait place à l'opérateur du laplacien
    MyMesh::Point laplacien_cot_v = 1/(2*aire) * somme;
    //laplacien_cot_v.normalize();

    return laplacien_cot_v;
}

MyMesh::Point MainWindow::calc_vector_v_vi(MyMesh* _mesh, VertexHandle v, VertexHandle vi){
    // On recupere les coordonnées de v et du voisin vi
    MyMesh::Point coords_v = _mesh->point(v);
    MyMesh::Point coords_vi = _mesh->point(vi);
    // on calcul le vecteur vvi
    MyMesh::Point vecteur = coords_vi - coords_v;

    return vecteur;
}
double MainWindow::calcul_poids_cot(MyMesh *_mesh, VertexHandle v, VertexHandle vi){
    // Un fonction de la classe MyMesh permet de calculer l'angle entre deux halfedges
    // Pour cela nous devons d'abords recupérer une des halfedge lié à l'edge vvi;
    HalfedgeHandle half_edge;
    for(MyMesh::VertexOHalfedgeIter he = _mesh->voh_iter(v); he.is_valid(); he++){
        // On s'assure que l'halfedge he pointe en direction de vi;
        if(_mesh->to_vertex_handle(he) == vi){
            half_edge = *he;
            break;
        }
    }

    // Calcul de alpha
    // Nous avons une halfedge qui correspond à l'edge vvi et a une des deux faces
    // Nous devons donc calculer l'angle avec les deux autres halfedges
    MyMesh::Scalar alpha = 0;
    HalfedgeHandle next = _mesh->next_halfedge_handle(half_edge);
    alpha = _mesh->calc_sector_angle(next);

    // Maintenant nous allons faire la meme chose avec l'halfedge opposé pour beta
    MyMesh::Scalar beta = 0;
    half_edge = _mesh->opposite_halfedge_handle(half_edge);
    next = _mesh->next_halfedge_handle(half_edge);
    beta = _mesh->calc_sector_angle(next);

    // On calcul ensuite la somme des cot de beta et alpha;
    MyMesh::Scalar somme = cot(alpha) + cot(beta);

    return somme;
}

MyMesh::Scalar MainWindow::neighboring_faces_area(MyMesh* _mesh, VertexHandle v){
    // On veut récuperer la somme du tiers des faces voisines au point v
    MyMesh::Scalar somme = 0;
    // On itere donc sur les faces voisines
    MyMesh::Scalar aire = 0;
    for(MyMesh::VertexFaceIter f = _mesh->vf_iter(v); f.is_valid(); f++){
        // Pour  calculer l'aire avec OpenMesh, il faut recuperer l'halfedge associé à la face
        HalfedgeHandle heh = _mesh->halfedge_handle(*f);
        // Puis calculer l'aire avec une fonction fourni par notre classe MyMesh
        aire = _mesh->calc_sector_area(heh);
        somme += aire;
    }

    // on divise ensuite la somme par trois
    somme = somme /3;
    return somme;
}


Matrix<MyMesh::Scalar, Dynamic, Dynamic> MainWindow::matrice_diag(MyMesh* _mesh)//Dynamic permet dene pas donner de taille fixe
{
    int nb_sommets = _mesh->n_vertices();
    Matrix<MyMesh::Scalar, Dynamic , Dynamic> d;
    d.resize(nb_sommets, nb_sommets);//Initialisation de la taille de la matrice
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end();v++)
    {
        VertexHandle vh = *v;
        MyMesh::Scalar aire = neighboring_faces_area(_mesh, vh);
        for(int i = 0; i < nb_sommets; i++)
        {
            for(int j =0; j < nb_sommets; j++)
            {
                if(i==j)
                    d(i,i) = 1/(2*(aire));
                else
                    d(i,j) = 0;
            }
        }
    }
    return d;
}//matrice_diag : renvoie la matrice diagonale du calcul de la matrice de Laplace Beltrami


Matrix<MyMesh::Scalar, Dynamic, Dynamic> MainWindow::matrice_adj(MyMesh* _mesh)//Dynamic permet de ne pas donner de taille fixe
{
    int nb_sommets = _mesh->n_vertices();
    Matrix<MyMesh::Scalar, Dynamic , Dynamic> m;
    m.resize(nb_sommets, nb_sommets);
    for(int i=0; i <nb_sommets; i++)
    {
        for(int j = 0; j<nb_sommets; j++)
        {
            m(i,j)=0;
        }
    }//Initialisation des valeurs de la matrice à 0
    MyMesh::Scalar poids;
    int cpt = 0;//Ce compteur sert a savoir ou nous sommes dans la matrice
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end();v++)
    {
        poids = 0;
        VertexHandle vh = *v;
        for(MyMesh::VertexVertexIter vi = _mesh->vv_iter(vh); vi.is_valid(); vi++)
        {
            VertexHandle vhi = *vi;
            poids += calcul_poids_cot(_mesh, vh, vhi);//calcul le poids par rapport a deux sommets adjacents et somme sur tous les sommets voisins de v
        }
        m(cpt,cpt) = -poids;
        cpt++;
    }
    return m;
}//matrice_adj : renvoie la matrice d'adjacence de Laplace Beltrami

Matrix<MyMesh::Scalar, Dynamic, Dynamic> MainWindow::matrice_Laplace_Beltrami(MyMesh* _mesh)
{
    int nb_sommets = _mesh->n_vertices();
    Matrix<MyMesh::Scalar, Dynamic , Dynamic> m;
    m.resize(nb_sommets, nb_sommets);//Initialisation de m
    m = matrice_adj(_mesh);//calcul de la matrice d'adjacence
    Matrix<MyMesh::Scalar, Dynamic , Dynamic> d;
    d.resize(nb_sommets, nb_sommets);//Initialisation de d
    d = matrice_diag(_mesh);//Calcul de la matrice diag
    return d*m;
}//matrice_Laplace_Beltrami : renvoie le produit de la matrice diag et de la matrice d'adjacence ce qui donne la matrice de Laplace Beltrami




void MainWindow::showSelections(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! ****
     * cette fonction utilise les vatiables de sélection vertexSelection, edgeSelection et faceSelection
     * qui sont les ID des élements sélectionnés et qui sont égales à -1 si la sélection est vide
     */

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}


void MainWindow::showSelectionsNeighborhood(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! ****
     * cette fonction utilise les vatiables de sélection vertexSelection, edgeSelection et faceSelection
     * qui sont les ID des élements sélectionnés et qui sont égales à -1 si la sélection est vide
     * et affiche en plus le voisinage de chaque sélection :
     *    - les faces voisines les faces
     *    - les faces adjacentes pour les arêtes
     *    - les arêtes incidentes pour les sommets
     */

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}



void MainWindow::showBorder(MyMesh* _mesh)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    /* **** à compléter ! **** */

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}


void MainWindow::showPath(MyMesh* _mesh, int v1, int v2)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // point de départ et point d'arrivée en vert et en gros
    _mesh->set_color(_mesh->vertex_handle(v1), MyMesh::Color(0, 255, 0));
    _mesh->set_color(_mesh->vertex_handle(v2), MyMesh::Color(0, 255, 0));
    _mesh->data(_mesh->vertex_handle(v1)).thickness = 12;
    _mesh->data(_mesh->vertex_handle(v2)).thickness = 12;

    /* **** à compléter ! **** */

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/* **** fin de la partie à compléter **** */


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);
    clone = mesh;
    // on affiche le maillage
    displayMesh(&mesh);
}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */

// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, DisplayMode mode)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(mode == DisplayMode::TemperatureMap)
    {
        QVector<float> values;
        for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
            values.append(fabs(_mesh->data(*curVert).value));
        qSort(values);

        float range = values.at(values.size()*0.8);

        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::Normal)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(mode == DisplayMode::ColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->data(*fvIt).faceShadingColor[0]; triCols[3*i+1] = _mesh->data(*fvIt).faceShadingColor[1]; triCols[3*i+2] = _mesh->data(*fvIt).faceShadingColor[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}



void MainWindow::on_operateur_clicked()
{
    operateur_laplace_beltrami(&clone, COTANGENTE, h ,_y);
    displayMesh(&clone);
}

void MainWindow::on_pushButton_clicked()
{
    displayMesh(&mesh);

}

void MainWindow::on_doubleSpinBox_valueChanged(double arg1)
{
    factor_change = true;
    h = -arg1;
}

void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
    factor_change = false;
    _y= arg1;
}

void MainWindow::on_scale_clicked()
{
    scale(&clone);
    scale(&mesh);
    displayMesh(&clone);
}

void MainWindow::on_Unscale_clicked()
{
    unscale(&clone);
    unscale(&mesh);
    displayMesh(&clone);
}

void MainWindow::on_pushButton_2_clicked()
{
    operateur_laplace_beltrami(&clone, UNIFORME, h ,_y);
    displayMesh(&clone);
}
