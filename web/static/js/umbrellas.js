$(function(){
    
    // a Parameter is just an id and a value
    window.Parameter = Backbone.Model.extend({ });
    
    // Each replica contains a ParameterSet
    window.ParameterSet = Backbone.Collection.extend({
        model: Parameter,
    });

    // Replica model
    window.Replica = Backbone.Model.extend({
        initialize: function(data) {
            this.set({id: data.name});
            this.parameters = new ParameterSet;
            this.parameters.url = '/replicas/' + this.id + '/parameters';
            // create all of the Parameter objects
            this.parameters.add(_.map(data.parameters, function(value, id){ return {id:id, value:value} }));
        },
        
        getSeries: function() {
            if (this.parameters.get("coordinate")) {
                return { name: this.id, data: [[this.parameters.get("coordinate").get("value"), 0.0]] };          
            }
            return false;
        },
    });

    // ReplicaSet
    window.ReplicaSet = Backbone.Collection.extend({
        model: Replica,
        url: '/replicas',
    });

    window.Replicas = new ReplicaSet;
    
    // Parameter view
    window.ParameterView = Backbone.View.extend({
        tagName: "div",
        className: "parameter",
        template: _.template($('#parameter-template').html()),
        
        events: {
            "click .edit-button": "edit",
            "click .save-button": "save",
            "click .cancel-button": "close",
        },
        
        initialize: function() {
            _.bindAll(this, 'render', 'edit', 'save', 'close');
            this.model.bind('change', this.render);
            this.model.view = this;
        },
        
        edit: function() {
            $(this.el).addClass("editing");
            this.input.focus();
        },
        
        save: function() {
            this.model.save({value: this.input.val()});
            this.close();
        },
        
        close: function() {
            this.setContent();
            $(this.el).removeClass("editing");
        },
        
        render: function() {
            $(this.el).html(this.template(this.model.toJSON()));
            this.setContent();
            return this;
        },
        
        setContent: function() {
            var content = this.model.get('value');
            this.$('.value').text(content);
            this.input = this.$('.value-input');
            this.input.val(content);
        },
        
    });
    
    // Replica view
    window.ReplicaView = Backbone.View.extend({
        tagName: "div",
        className: "replica",
        template: _.template($('#replica-template').html()),
        
        events: {
        },
        
        initialize: function() {
            _.bindAll(this, 'render');
            // Bind the view to the model (one-to-one)
            this.model.bind('change', this.render);
            this.model.view = this;              
        },
        
        // Render the contents of a single replica.
        render: function() {
            $(this.el).html(this.template(this.model.toJSON()));
            // add the parameter sub-templates into this one
            var self = this;
            _.each(self.model.parameters.models, function(p) { self.$('.parameter-list').append(new ParameterView({model: p}).render().el); });
            return this;
        },
        
    });
    
    // The Application

    // Our overall **AppView** is the top-level piece of UI.
    window.AppView = Backbone.View.extend({
        el: $("#umbrellas-app"),

        events: {
        },

        initialize: function() {
            _.bindAll(this, 'addOne', 'addAll', 'render');
            Replicas.bind('add',     this.addOne);
            Replicas.bind('refresh', this.addAll);
            Replicas.bind('all',     this.render);
        },

        render: function() {
        },
        
        // Add a single replica to the list by creating a view for it, and
        // appending its element to the `<ul>`.
        addOne: function(replica) {
            var view = new ReplicaView({model: replica});
            this.$("#replica-list").append(view.render().el);
            series = replica.getSeries()
            if (series) {
                this.chart.addSeries(series);
            }
        },
        
        addAll: function() {
            Replicas.each(this.addOne);
        }

    });

    // Finally, we kick things off by creating the **App**.
    window.App = new AppView;

});