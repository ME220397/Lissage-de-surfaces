#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<cmath>
#include<QVector3D>
#include "materiel_courbures/courbures.h"
/* **** début de la partie à compléter **** */
Courbures *courb;
double cot(double teta){
    return std::cos(teta)/std::sin(teta);
}

void scale(MyMesh * _mesh){
    Vec3f new_coords;
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        new_coords = _mesh->point(*v);
        _mesh->set_point(*v, new_coords*24);
    }
}

void unscale(MyMesh * _mesh){
    Vec3f new_coords;
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        new_coords = _mesh->point(*v);
        _mesh->set_point(*v, new_coords*1/24);
    }
}


void MainWindow::flou_de_diffusion(MyMesh *_mesh, VertexHandle v, MyMesh::Point fv, double f){
    // On calcul le nouveau point v
    MyMesh::Point new_v = _mesh->point(v) + f * fv;
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
    }else{
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
    double n_voisin = 0;
    for(MyMesh::VertexVertexIter vi = _mesh->vv_iter(v); vi.is_valid(); vi++){
        // 3. Calcul des du vecteur directeur vvi
        direction = calc_vector_v_vi(_mesh, v, *vi);
        u = direction;
        //u.normalize();
        somme += u;
        n_voisin++;
    }
    MyMesh::Point laplacien_uni_v = (1.0/n_voisin) * somme;
    //laplacien_cot_v.normalize();

    return laplacien_uni_v;

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
    MyMesh::Point laplacien_cot_v = 1/(2.0*aire) * somme;
    qDebug() << laplacien_cot_v[0] << laplacien_cot_v[1] << laplacien_cot_v[2];
    laplacien_cot_v.normalize();
    qDebug() << laplacien_cot_v[0] << laplacien_cot_v[1] << laplacien_cot_v[2];

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
double norm(Vec3f v){
    return sqrt(pow(v[0], 2) + pow(v[0], 2) +pow(v[0], 2));
}
double produit_scalaire(Vec3f u, Vec3f v){
    double ux = u[0], uy = u[1], uz = u[2];
    double vx = v[0], vy = v[1], vz = v[2];

    return ux*vx + uy*vy + uz*vz;
}
double calc_angle(MyMesh::Point v, MyMesh::Point vi, MyMesh::Point p){
    double vx = v[0], vy = v[1], vz = v[2];
    double vix = vi[0], viy = vi[1], viz = vi[2];
    double px = p[0], py = p[1], pz = p[2];

    QVector3D u(vx - px, vy - py, vz - pz);
    QVector3D w(vix - px, viy - py, viz - pz);
    u.normalize();
    w.normalize();
    double cos = QVector3D::dotProduct(u, w);
    double angle = acos(cos);
    return angle;

}
double MainWindow::calcul_poids_cot(MyMesh *_mesh, VertexHandle v, VertexHandle vi){
    // Un fonction de la classe MyMesh permet de calculer l'angle entre deux halfedges
    // Pour cela nous devons d'abords recupérer une des halfedge lié à l'edge vvi;
    HalfedgeHandle half_edge;
    for(MyMesh::VertexOHalfedgeIter he = _mesh->voh_iter(v); he.is_valid(); he++){
        // On s'assure que l'halfedge he pointe en direction de vi;
        if(_mesh->to_vertex_handle(*he) == vi){
            half_edge = *he;
            break;
        }
    }
    // Calcul de alpha
    // Nous avons une halfedge qui correspond à l'edge vvi et a une des deux faces
    // Nous devons donc calculer l'angle avec les deux autres halfedges
    MyMesh::Scalar alpha = 0;
    HalfedgeHandle next = _mesh->next_halfedge_handle(half_edge);
    VertexHandle vh = _mesh->to_vertex_handle(next);
    MyMesh::Point right = _mesh->point(vh);
    alpha = calc_angle(_mesh->point(v), _mesh->point(vi), right);

    // Maintenant nous allons faire la meme chose avec l'halfedge opposé pour beta
    MyMesh::Scalar beta = 0;
    half_edge = _mesh->opposite_halfedge_handle(half_edge);
    next = _mesh->next_halfedge_handle(half_edge);
    vh = _mesh->to_vertex_handle(next);
    MyMesh::Point left = _mesh->point(vh);
    beta = calc_angle(_mesh->point(v), _mesh->point(vi), left);

    // On calcul ensuite la somme des cot de beta et alpha;
    MyMesh::Scalar somme = cot(alpha) + cot(beta);

    return somme;
}

double calcul_aire(MyMesh::Point p[]){
    MyMesh::Point v1 = p[0];
    MyMesh::Point v2 = p[1];
    MyMesh::Point v3 = p[2];

    // on recupere le vecteur v1v2
    MyMesh::Point u = v2 - v1;
    // on recupere le vecteur v1v3
    MyMesh::Point v = v3 - v1;

    double ux = u[0], uy = u[1], uz = u[2];
    double vx = v[0], vy = v[1], vz = v[2];

    // on calcul les determinant du produit vectoriel
    volatile double i = (uy*vz - vy*uz);
    volatile double j = (ux*vz - vx*uz);
    volatile double k = (ux*vy - vx*uy);

    volatile double area2 = i - j + k;
    volatile double area = abs(area2)/2.0;
    return area;

}
MyMesh::Scalar MainWindow::neighboring_faces_area(MyMesh* _mesh, VertexHandle v){
    // On veut récuperer la somme du tiers des faces voisines au point v
    MyMesh::Scalar somme = 0;
    // On itere donc sur les faces voisines
    MyMesh::Scalar aire = 0;
    MyMesh::Point points[3];
    int cpt= 0;
    for(MyMesh::VertexFaceIter f = _mesh->vf_iter(v); f.is_valid(); f++){
        cpt = 0;
        // Pour  calculer l'aire avec OpenMesh, il faut recuperer l'halfedge associé à la face
        for(MyMesh::FaceVertexIter fv_it = _mesh->fv_iter(*f); fv_it.is_valid(); fv_it++){
            points[cpt++] = _mesh->point(*fv_it);
        }
        // Puis calculer l'aire avec une fonction fourni par notre classe MyMesh
        aire = calcul_aire(points);
        somme += aire;
    }

    // on divise ensuite la somme par trois
    return somme/3;
}

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

    if (mode == DisplayMode::VertexColorShading)
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fvIt)[0]; triCols[3*i+1] = _mesh->color(*fvIt)[1]; triCols[3*i+2] = _mesh->color(*fvIt)[2];
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

void MainWindow::calc_courbure(MyMesh *_mesh){
    _mesh->update_normals();
    _mesh->request_vertex_colors() ;

    courb = new Courbures(*_mesh) ;
    resetAllColorsAndThickness(_mesh);
    courb->compute_KH();
    courb->set_K_colors(H);
    displayMesh(&mesh,DisplayMode::ColorShading);
}

void set_white_color(MyMesh * _mesh){
    for(MyMesh::VertexIter v = _mesh->vertices_begin(); v != _mesh->vertices_end(); v++){
        _mesh->set_color(*v, MyMesh::Color(150,150,150));
    }
}
void MainWindow::on_Unscale_clicked()
{
    calc_courbure(&clone);
    displayMesh(&clone, DisplayMode::VertexColorShading);
}

void MainWindow::on_pushButton_2_clicked()
{
    operateur_laplace_beltrami(&clone, UNIFORME, h ,_y);
    displayMesh(&clone);
}
