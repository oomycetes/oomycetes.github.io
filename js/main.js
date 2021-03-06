var species_dict = {};
species_dict["ALCA"] = "Albugo candida";
species_dict["ALLA"] = "Albugo laibachii";
species_dict["APAS"] = "Aphanomyces astaci";
species_dict["APIN"] = "Aphanomyces invadans";
species_dict["HYAP"] = "Hyaloperonospora arabidopsis";
species_dict["PHAG"] = "Phytophthora agathidicida";
species_dict["PHCA"] = "Phytophthora capsici";
species_dict["PHCI"] = "Phytophthora cinnamomi";
species_dict["PHCR"] = "Phytophthora cryptogea";
species_dict["PHFR"] = "Phytophthora fragariae";
species_dict["PHIF"] = "Phytophthora infestans";
species_dict["PHKE"] = "Phytophthora kernoviae";
species_dict["PHLA"] = "Phytophthora lateralis";
species_dict["PHMU"] = "Phytophthora multivora";
species_dict["PHNI"] = "Phytophthora nicotianae";
species_dict["PHPA"] = "Phytophthora parasitica";
species_dict["PHPI"] = "Phytophthora pinifolia";
species_dict["PHPL"] = "Phytophthora pluvialis";
species_dict["PHPS"] = "Phytophthora pisi";
species_dict["PHRA"] = "Phytophthora ramorum";
species_dict["PHRU"] = "Phytophthora rubi";
species_dict["PHSO"] = "Phytophthora sojae";
species_dict["PHTT"] = "Phytophthora taxon totara";
species_dict["PIAP"] = "Pilasporangium apinafurcum";
species_dict["PLHA"] = "Plasmopara halstedii";
species_dict["PLVI"] = "Plasmopara viticola";
species_dict["PYAP"] = "Pythium aphanidermatum";
species_dict["PYAR"] = "Pythium arrhenomanes";
species_dict["PYIN"] = "Pythium insidiosum";
species_dict["PYIR"] = "Pythium irregulare";
species_dict["PYIW"] = "Pythium iwayami";
species_dict["PYOL"] = "Pythium oligandrum ";
species_dict["PYUS"] = "Pythium ultimum var. sporangiiferum";
species_dict["PYUU"] = "Pythium ultimum var. ultimum";
species_dict["PYVX"] = "Phytopythium vexans";
species_dict["SADI"] = "Saprolegnia diclina";
species_dict["SAPA"] = "Saprolegnia parasitica";

var order_dict = {};
order_dict["ALCA"] = "Albuginales";
order_dict["ALLA"] = "Albuginales";
order_dict["APAS"] = "Saprolegniales";
order_dict["APIN"] = "Saprolegniales";
order_dict["HYAP"] = "Peronosporales";
order_dict["PHAG"] = "Peronosporales";
order_dict["PHCA"] = "Peronosporales";
order_dict["PHCI"] = "Peronosporales";
order_dict["PHCR"] = "Peronosporales";
order_dict["PHFR"] = "Peronosporales";
order_dict["PHIF"] = "Peronosporales";
order_dict["PHKE"] = "Peronosporales";
order_dict["PHLA"] = "Peronosporales";
order_dict["PHMU"] = "Peronosporales";
order_dict["PHNI"] = "Peronosporales";
order_dict["PHPA"] = "Peronosporales";
order_dict["PHPI"] = "Peronosporales";
order_dict["PHPL"] = "Peronosporales";
order_dict["PHPS"] = "Peronosporales";
order_dict["PHRA"] = "Peronosporales";
order_dict["PHRU"] = "Peronosporales";
order_dict["PHSO"] = "Peronosporales";
order_dict["PHTT"] = "Peronosporales";
order_dict["PIAP"] = "Pythiales";
order_dict["PLHA"] = "Peronosporales";
order_dict["PLVI"] = "Peronosporales";
order_dict["PYAP"] = "Pythiales";
order_dict["PYAR"] = "Pythiales";
order_dict["PYIN"] = "Pythiales";
order_dict["PYIR"] = "Pythiales";
order_dict["PYIW"] = "Pythiales";
order_dict["PYOL"] = "Pythiales";
order_dict["PYUS"] = "Pythiales";
order_dict["PYUU"] = "Pythiales";
order_dict["PYVX"] = "Peronosporales";
order_dict["SADI"] = "Saprolegniales";
order_dict["SAPA"] = "Saprolegniales";

var current_id;
var current_sequence;

function nodeClicked(node) {
	current_id = node.data.node.id;
	$('#modal-label').text(node.data.node.id);
	$('#protein-id').text(node.data.node.id);
	$('#species').text(species_dict[node.data.node.id.substr(0,4)]);
	$('#genus').text(species_dict[node.data.node.id.substr(0,4)].split(" ")[0]);
	$('#order').text(order_dict[node.data.node.id.substr(0,4)]);
	$('#degree').text(parseInt(graph.graph.degree(node.data.node.label)));

	// Get info regarding methods used to identify candidate RxLR
	if(window.location.href.includes("rxlrs")) {
		if(node.data.node.attributes.Win == 1) {
			$('#win-method').text("Yes")
		}
		else {
			$('#win-method').text("No")	
		}
		if(node.data.node.attributes.Regex == 1) {
			$('#regex-method').text("Yes")
		}
		else {
			$('#regex-method').text("No")	
		}
		if(node.data.node.attributes.HMM == 1) {
			$('#hmm-method').text("Yes")
		}
		else {
			$('#hmm-method').text("No")	
		}
		if(node.data.node.attributes.BLAST == 1) {
			$('#blast-method').text("Yes")
		}
		else {
			$('#blast-method').text("No")	
		}
		if(node.data.node.attributes.WYLdomain == 'Y') {
			$('#wyl-domain').text("Yes")
		}
		else {
			$('#wyl-domain').text("No")	
		}
	}


	if(window.location.href.includes("chitinase")) {
		getSequenceChitinases(current_id);
	} else if(window.location.href.includes("nlps")) {
		getSequenceNLPs(current_id);
	} else if(window.location.href.includes("rxlrs")) {
		getSequenceRxLRs(current_id);
	} else {
		console.log("Error!");
	}

	$('#modal').modal('show');
}

function getSequenceChitinases(id) {
	file_path = 'data/fasta/chitinases/' + id + '.fasta';
	getSequence(id, file_path) 
}

function getSequenceNLPs(id) {
	file_path = 'data/fasta/nlps/' + id + '.fasta';
	getSequence(id, file_path) 
}

function getSequenceRxLRs(id) {
	file_path = 'data/fasta/rxlrs/' + id + '.fasta';
	getSequence(id, file_path) 
}

function getSequence(id, file_path) {
	$.get(file_path, function(data) {
	    current_sequence = data.split("\n").slice(1).join("");
	    $('#protein-length').text(current_sequence.length);
	    $('#sequence').text(current_sequence);
	}, 'text');
};

$("#button-close_help-modal1").click(function() {
	$("#modal-help").modal("hide");
})

$("#button-close_help-modal2").click(function() {
	$("#modal-help").modal("hide");
})

$("#blast").click(function() {
	window.open("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY=" + current_sequence)
});

// Filter the network, only showing sequences that have ID/species/order matching filter text box
$('#filter').on('input', function() {
	filterNetwork();
});

function filterNetwork() {
	var $this = $(this);
    var delay = 500; // Delay for 0.5 seconds

    var n_visible = graph.graph.nodes().length // Use to count how many nodes are visible

    // Check if this is the RxLR page, if so need to consider additional filtering options
    if (typeof rxlr_page !== 'undefined') {
		// This is the RxLR page
		clearTimeout($this.data('timer'));
	    $this.data('timer', setTimeout(function(){
	        $this.removeData('timer');

	        filter_value = $('#filter').val().trim().toUpperCase();

			if( $("#win-check").is(":checked")) {
				win_check = true;
			}
			else {
				win_check = false;
			}

			if( $("#regex-check").is(":checked")) {
				regex_check = true;
			}
			else {
				regex_check = false;
			}

			if( $("#hmm-check").is(":checked")) {
				hmm_check = true;
			}
			else {
				hmm_check = false;
			}

			if( $("#blast-check").is(":checked")) {
				blast_check = true;
			}
			else {
				blast_check = false;
			}

			if( $("#wyl-check").is(":checked")) {
				wyl_check = true;
			}
			else {
				wyl_check = false;
			}

			// Get nodes and hide those that don't match filter or method checkboxes
			nodes = graph.graph.nodes();

			// If filter textbox is empty show all nodes that match methods option
			if(filter_value.length == 0) {
				for(i = 0; i < nodes.length; i++) {
					methods = getNodeValues(nodes[i]);
					if( win_check && methods['win'] == true ) {
						nodes[i].hidden = false;
					}
					else if ( regex_check && methods['regex'] == true ) {
						nodes[i].hidden = false;
					}
					else if ( hmm_check && methods['hmm'] == true ) {
						nodes[i].hidden = false;
					}
					else if ( blast_check && methods['blast'] == true ) {
						nodes[i].hidden = false;
					}
					else if ( wyl_check && methods['wyl'] == true ) {
						nodes[i].hidden = false;
					}
					else {
						nodes[i].hidden = true;
						n_visible -= 1;
					}
				}
				$("#visible_nodes").text(n_visible);
			}
			else {
				for(i = 0; i < nodes.length; i++) {
					// Find out if nodes hits at least on of the selected/checked methods
					methods = getNodeValues(nodes[i]);
					match_method = false;

					if( win_check && methods['win'] == true ) {
						match_method = true;
					}
					else if ( regex_check && methods['regex'] == true ) {
						match_method = true;
					}
					else if ( hmm_check && methods['hmm'] == true ) {
						match_method = true;
					}
					else if ( blast_check && methods['blast'] == true ) {
						match_method = true;
					}
					else if ( wyl_check && methods['wyl'] == true ) {
						match_method = true;
					}

					if(nodes[i].id.includes(filter_value) || idToSpecies(nodes[i].id).includes(filter_value) || idToOrder(nodes[i].id).includes(filter_value)) {
						if (match_method) {
							nodes[i].hidden = false;
						}
						else {
							nodes[i].hidden = true;
							n_visible -= 1;
						}
					}
					else {
						nodes[i].hidden = true;
						n_visible -= 1
					}
				}
				$("#visible_nodes").text(n_visible);
			}
			graph.refresh({ skipIndexation: true});
	    }, delay));

	}
	else {
		clearTimeout($this.data('timer'));
	    $this.data('timer', setTimeout(function(){
	        $this.removeData('timer');

	        filter_value = $('#filter').val().trim().toUpperCase();
			// Get nodes and hide those that don't match filter
			nodes = graph.graph.nodes();

			// If filter textbox is empty show all nodes
			if(filter_value.length == 0) {
				for(i = 0; i < nodes.length; i++) {
					nodes[i].hidden = false;
				}
			}
			else {
				for(i = 0; i < nodes.length; i++) {
					if(nodes[i].id.includes(filter_value) || idToSpecies(nodes[i].id).includes(filter_value) || idToOrder(nodes[i].id).includes(filter_value)) {
						nodes[i].hidden = false;
					}
					else {
						nodes[i].hidden = true;
						n_visible -= 1;
					}
				}
			}
			$("#visible_nodes").text(n_visible);
			graph.refresh({ skipIndexation: true});
	    }, delay));
	}
}

function idToSpecies(input_id) {
	return(species_dict[input_id.substr(0,4)]).toUpperCase();
}

function idToOrder(input_id) {
	return(order_dict[input_id.substr(0,4)]).toUpperCase();
}

// Returns the methods/features that identified the node as a putative RxLR
// Returns boolean array of Win, Regex, HMM, BLAST and WYL
function getNodeValues(node) {
	result = [];
	if(node.attributes.Win == "1") {
		result["win"] = true;

	}
	else {
		result["win"] = false;
	}
	if(node.attributes.Regex == "1") {
		result['regex'] = true;
	}
	else {
		result['regex'] = false;
	}
	if(node.attributes.HMM == "1") {
		result['hmm'] = true;
	}
	else {
		result['hmm'] = false;
	}
	if(node.attributes.BLAST == "1") {
		result['blast'] = true;
	}
	else {
		result['blast'] = false;
	}
	if(node.attributes.WYLdomain == "Y") {
		result['wyl'] = true;
	}
	else {
		result['wyl'] = false;
	}
	return result;
}

$("#win-check").click(function() {
	filterNetwork();
})

$("#regex-check").click(function() {
	filterNetwork();
})

$("#hmm-check").click(function() {
	filterNetwork();
})

$("#blast-check").click(function() {
	filterNetwork();
})

$("#wyl-check").click(function() {
	filterNetwork();
})
