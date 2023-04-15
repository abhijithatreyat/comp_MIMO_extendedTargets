function [outputArg1,outputArg2] = surface_from_scatter(ptCloud,disp_ptCloud)
    %% Making Surface Plots From Scatter Data
    % How do you turn a collection of XYZ triplets into a surface plot? This is
    % the most frequently asked 3D plotting question that I got when I was in
    % Tech Support.
    %% Load the data
    
    %%
    if(disp_ptCloud)
        figure;
        plot3(x,y,z,'.-');
        title("Input to surface_from_scatter");
    end
    %% Little triangles
    % The solution is to use Delaunay triangulation. Let's look at some
    % info about the "tri" variable.
    tri = delaunay(x,y);
    plot(x,y,'.')
    %%
    % How many triangles are there?
    [r,c] = size(tri);
    disp(r)
    %% Plot it with TRISURF
    if(disp_ptCloud)
        figure;
        h = trisurf(tri, x, y, z);
        title("Trisurf output");
        axis vis3d;
        axis off;
    end
    %% Clean it up
    
    l = light('Position',[-50 -15 29]);
    set(gca,'CameraPosition',[208 -50 7687]);
    lighting phong
    shading interp
    colorbar EastOutside

end

