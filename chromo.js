/*eslint-disable no-unused-vars */
/*eslint-env es6 */
/*global document */

//
// Author: H. del Risco
//

class Chromo {
    constructor(name, bases, maxBases) {
        this.name = name;
        this.bases = bases;
        this.maxBases = maxBases;
    }

    // caller must decide where to place the svg as to have the proper margins in the overall display
    createElement(width, height, bands) {
        var ppb = width / this.maxBases,
            n = bands.length,
            preLast;

        // create SVG element and add bands
        var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("shape-rendering", "geometricPrecision");
        svg.setAttribute("width", "".concat(width));
        svg.setAttribute("height", "".concat(height));
        var element, title, //acenEnd,
            acen = -1;
        var html = "",
            strfill = "none",
            endfill = "none",
            i;
        for(i = 0; i < n; i++) {
            // add band
            var band = bands[i];            
            var x = band.strpos * ppb;
            element = document.createElementNS("http://www.w3.org/2000/svg", "rect");
            var wsub = 0;
            if(i === 0) {
                x += 2;
                wsub = 2;
                strfill = band.color;
            }
            else if(i === (n - 1)) {
                wsub = 2;
                endfill = band.color;
            }
            var size = Math.abs(band.endpos - band.strpos);
            var bw = size * ppb - wsub;
            if(band.stain === "acen" && acen === (i-1)) {
                x -= 1;
                bw += 1;
            }
            if(bw < 0)
                continue;
            element.setAttribute("x", "".concat(x));
            element.setAttribute("y", "0");
            element.setAttribute("width", "".concat(bw));
            element.setAttribute("height", "".concat(height));
            element.style.fill = band.color;
            title = document.createElementNS("http://www.w3.org/2000/svg", "title");
            title.textContent = "Chromosome: " + this.name + "\nBand: " + band.name + "\nStain: " + band.stain + "\nPosition: " + band.strpos + " to " + band.endpos;
            element.appendChild(title);
            svg.appendChild(element);
            if(band.stain === "acen") {
                if(acen === -1) {
                    acen = i;
                }
                else {
                    if(acen === (i-1)) {
                        //addLine = false;
                        element = document.createElementNS("http://www.w3.org/2000/svg", "line");
                        element.setAttribute("x1", "".concat(x));
                        element.setAttribute("y1", "0"); // should offset chr to 1 to allow using 0 here to get rid of the stroke
                        element.setAttribute("x2", "".concat(x));
                        element.setAttribute("y2", "".concat(height/2 - 5));
                        element.style.stroke = "white";
                        element.style.strokeWidth = "2.0";
                        svg.appendChild(element);

                        element = document.createElementNS("http://www.w3.org/2000/svg", "line");
                        element.setAttribute("x1", "".concat(x));
                        element.setAttribute("y1", "".concat(height/2 + 5));
                        element.setAttribute("x2", "".concat(x));
                        element.setAttribute("y2", "".concat(height));
                        element.style.stroke = "white";
                        element.style.strokeWidth = "2.0";
                        svg.appendChild(element);
                    }
                }
            }
        }
        
        width = Math.max(this.bases * ppb, 4);
        
        var path = "M2 0 C 0 0, 0 " + height + ", 2 " + height;
        element = document.createElementNS("http://www.w3.org/2000/svg", "path");
        element.setAttribute("d", path);
        element.style.stroke = "black";
        element.style.fill = strfill;            
        element.style.strokeWidth = "0.25";
        svg.appendChild(element);

        if(width > 4) {
            path = " M2 0 L " + (width - 3) + " 0";
            path += " M2 " + height + " L " + (width - 3) + " " + height;
            element = document.createElementNS("http://www.w3.org/2000/svg", "path");
            element.setAttribute("d", path);
            element.style.stroke = "black";
            element.style.fill = "none";            
            element.style.strokeWidth = "0.5";
            // check if no bands provided
            if(bands.length === 0) {
                // setup tooltip for the whole chromosome w/o bands - must fill to get it
                title = document.createElementNS("http://www.w3.org/2000/svg", "title");
                title.textContent = "Chromosome: " + this.name + "\nPosition: 0 to " + this.bases;
                element.style.fill = "white";
                element.appendChild(title);            
            }
            svg.appendChild(element);
        }
        
        path = "M" + (width - 3) + " 0 C " + (width - 1) + " 0, " + (width - 1) + " " + height + ", " + (width - 3) + " " + height;
        element = document.createElementNS("http://www.w3.org/2000/svg", "path");
        element.setAttribute("d", path);
        element.style.stroke = "black";
        element.style.fill = endfill;            
        element.style.strokeWidth = "0.25";
        svg.appendChild(element);
        return svg;
    }
    getSVG(width, height, bands) {
        var ppb = width / this.maxBases,
            n = bands.length,
            svg = "",
            preLast,
            element, title,
            acen = -1,
            html = "",
            strfill = "none",
            endfill = "none",
            i;
        
        // add svg element tag
        svg += '<svg width="' + width + '" height="' + height + '"';
        svg += ' shape-rendering="geometricPrecision">';

        // add tags for all bands
        for(i = 0; i < n; i++) {
            var band = bands[i];

            // add band
            var x = band.strpos * ppb;
            var wsub = 0;
            if(i === 0) {
                x += 2;
                wsub = 2;
                strfill = band.color;
            }
            else if(i === (n - 1)) {
                wsub = 2;
                endfill = band.color;
            }
            var size = Math.abs(band.endpos - band.strpos);
            var bw = size * ppb - wsub;
            if(band.stain === "acen" && acen === (i-1)) {
                x -= 1;
                bw += 1;
            }
            if(bw < 0)
                continue;
            svg += '<rect x="' + x + '" y="0" width="' + bw + '" height="' + height + '" style="fill:' + band.color + ';">';
            svg += '<title>Chromosome: ' + this.name + '\nBand: ' + band.name + '\nStain: ' + band.stain + '\nPosition: ' + band.strpos + ' to ' + band.endpos + '</title></rect>';
            
            if(band.stain === "acen") {
                if(acen === -1)
                    acen = i;
                else {
                    if(acen === (i - 1)) {
                        svg += '<line x1="' + x + '" y1="0" x2="' + x + '" y2="' + (height/2 - 5) + '" style="stroke:white; stroke-width:2;"/>';
                        svg += '<line x1="' + x + '" y1="' + (height/2 + 5) + '" x2="' + x + '" y2="' + height + '" style="stroke:white; stroke-width:2;"/>';
                    }
                }
            }
        }
        
        width = Math.max(this.bases * ppb, 4);
        svg += '<path d="M2 0 C 0 0, 0 ' + height + ', 2 ' + height + '" style="stroke:black; stroke-width:0.25; fill:' + strfill + ';"/>';

        if(width > 4) {
            svg += '<path d="M2 0 L ' + (width - 3) + ' 0 M2 ' + height + ' L ' + (width - 3) + ' ' + height + '" style="stroke:black; stroke-width:0.5; fill:none;"';
            // check if no bands provided
            if(bands.length === 0) {
                // setup tooltip for the whole chromosome w/o bands - must use fill to see it
                svg += '><title style="fill:white;">Chromosome: ' + this.name + '\nPosition: 0 to ' + this.bases + '</title></path>';
            }
            else
                svg += "/>"
        }
        
        svg += '<path d="M' + (width - 3) + ' 0 C ' + (width - 1) + ' 0, ' + (width - 1) + ' ' + height + ', ' + (width - 3) + ' ' + height + '" style="stroke:black; stroke-width:0.25; fill:' + endfill + ';"/>';
        svg += "</svg>"
        return svg;
    }
}

